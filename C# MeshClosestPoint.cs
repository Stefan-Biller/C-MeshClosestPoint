// Grasshopper Script Instance
#region Usings
using System;
using System.Linq;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Collections.Concurrent;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
#endregion

struct Timings
    {
        public double setup;
        public double execution;
        public double total;
    };

public class Script_Instance : GH_ScriptInstance
{
    #region Notes
    /* 
      Members:
        RhinoDoc RhinoDocument
        GH_Document GrasshopperDocument
        IGH_Component Component
        int Iteration

      Methods (Virtual & overridable):
        Print(string text)
        Print(string format, params object[] args)
        Reflect(object obj)
        Reflect(object obj, string method_name)
    */
    #endregion

    private void RunScript(
		Mesh M,
		List<Point3d> Pts,
		int layers,
		int nodeSize,
		int iterations,
		ref object P,
		ref object T,
		ref object t)
    {
        if (iterations < 1) return;

        layers = Math.Max(1, layers);
        nodeSize = Math.Max(1, nodeSize);
        Stopwatch sw = new Stopwatch();
        
        int PARALLEL_BATCH_SIZE = 32;
        var totalTimeMs = new Timings();

        IEnumerable<GH_Point> lastClosestPoints = Enumerable.Empty<GH_Point>();
        BVHTree lastTree = null;

        for (int iter = 0; iter < iterations; iter++)
        {
            // setup
            sw.Start();
            var closestPoints = new GH_Point[Pts.Count];
            var tree = new BVHTree(M, nodeSize, layers);
            var sortedPoints = new Point3d[Pts.Count];
            var sortedPointsIndexMap = new int[Pts.Count];
            TreeHelper.SortPointsByMorton(Pts, out sortedPoints, out sortedPointsIndexMap);
            sw.Stop();
            double setupMs = sw.ElapsedMilliseconds;
            totalTimeMs.setup += setupMs;
            sw.Reset();

            // execution
            sw.Start();
            var batches = Partitioner.Create(0, Pts.Count, PARALLEL_BATCH_SIZE);
            Parallel.ForEach(batches, (batch) => {
                int start = batch.Item1;
                int end = batch.Item2;

                double maxSearchDistance = double.PositiveInfinity;
                var lastClosestPoint = Point3d.Unset; // keep last cp ref to limit next search radius
                
                for (int i = start; i < end; i++)
                {
                    var pt = sortedPoints[i];
                    maxSearchDistance = lastClosestPoint.IsValid ? pt.DistanceTo(lastClosestPoint) : double.PositiveInfinity;
                    lastClosestPoint = tree.closestPoint(pt, maxSearchDistance);
                    closestPoints[sortedPointsIndexMap[i]] = new GH_Point(lastClosestPoint);
                }
                
            });
            sw.Stop();
            lastClosestPoints = closestPoints;
            double executionMs = sw.ElapsedMilliseconds;
            totalTimeMs.execution += executionMs;
            sw.Reset();
            

            double iterationTotalMs = setupMs + executionMs;
            Print($"Iteration: {iterationTotalMs} ms");
            t = iterationTotalMs;
            totalTimeMs.total += iterationTotalMs;
            lastTree = tree;
        }

        // output
        P = lastClosestPoints;
        T = lastTree.bbTree;

        Print("______total Time________");
        Print($"Setup: {totalTimeMs.setup} ms");
        Print($"Execution: {totalTimeMs.execution} ms");
        Print("______last result_______");
        Print($"Tree Level Checks: {string.Join(" | ", lastTree.layerChecks)}");
        Print($"Triangle Checks: {lastTree.triangleChecks.ToString()}");
        Print($"Triangle Checks / Point: {(lastTree.triangleChecks / Pts.Count).ToString()}");
        Print($"Total: {totalTimeMs.total} ms");
        Print("______Avg______________");
        Print($"Average: {totalTimeMs.total / iterations} ms");
    }
}

public class BVHTree
{
    // continuous triangle list for improved caching during evaluation
    public Triangle3d[] triangles;
    public Point3d[] triangleVertices;
    // bounding box tree - layers are ordered [coarsest, ...,  finest]
    public BoundingBox[][] bbTree;

    public readonly int chunkSize;
    private readonly int treeLayers;
    private readonly double tolerance;

    // debugging and performance eval
    public int[] layerChecks;
    public int triangleChecks;

    private const string OUT_OF_BOUNDS_ERROR = "tree depth out of bounds";

    public BVHTree (
        Mesh mesh, 
        int chunkSize, 
        int treeLayers,
        double tolerance = 1e-9 // Rhino mm tolerance
        ) 
    {
        treeLayers = treeLayers < 1 ? 1 : treeLayers;

        this.chunkSize = chunkSize < 1 ? 1 : chunkSize;
        this.tolerance = tolerance < 0 ? 0 : tolerance;
        this.treeLayers = treeLayers < 1 ? 1 : treeLayers;
        this.bbTree = new BoundingBox[treeLayers][];
        this.layerChecks = new int[treeLayers];

        // build continuous triangle list from mesh
        bool triangulateResult = mesh.Faces.ConvertQuadsToTriangles();
        if (!triangulateResult) throw new ArgumentException("invalid mesh, could not triangulate");
        this.triangles = new Triangle3d[mesh.Faces.Count];

        for (int i = 0; i < mesh.Faces.Count; i++)
        {
            MeshFace face = mesh.Faces.GetFace(i);
            this.triangles[i] = new Triangle3d(
                mesh.Vertices[face.A],
                mesh.Vertices[face.B],
                mesh.Vertices[face.C]
            );
        }
            
        // reorded triangles to reduce bounding boxes / improve locality
        // use Morton Codes to reorder triangles
        TreeHelper.SortTrianglesByMorton(this.triangles);

        this.triangleVertices = this.triangles.SelectMany(t => new [] { t.A, t.B, t.C }).ToArray();


        // build bb tree
        // root layer - bb encapsulating triangle chunks
        var leafCount = (int) Math.Ceiling((double) this.triangles.Length / chunkSize);
        this.bbTree[this.treeLayers - 1] = new BoundingBox[leafCount];
        for (int i = 0; i < leafCount; i++)
        {
            this.bbTree[this.treeLayers - 1][i] = BVHTree.BoundingBoxFromTriangles(
                        this.triangles, 
                        i * chunkSize, 
                        chunkSize
                    );
        }
        // upper layers -- bb encapsulating bb chunks
        for (int d = this.treeLayers - 2; d >= 0; d--)
        {
            var lowerLayer = this.bbTree[d + 1];
            int layerCount = (int) Math.Ceiling((double) lowerLayer.Length / chunkSize);
            this.bbTree[d] = new BoundingBox[layerCount];
            for (int i = 0; i < layerCount; i++)
            {
                this.bbTree[d][i] = BVHTree.BoundingBoxFromBoundingBoxes(
                    lowerLayer,
                    i * chunkSize, 
                        chunkSize
                );
            }
            
        }

    }


    public Point3d closestPoint(Point3d testPoint, double maxDistance = double.PositiveInfinity) 
    {
        Point3d closestPoint = Point3d.Unset;
        double distanceSqrd = double.IsInfinity(maxDistance) ? maxDistance : maxDistance * maxDistance;

        if (double.IsInfinity(distanceSqrd))
        {
            // pick an initial distance from some samples in the tree
            this.closestPointAtTriangleChunk(0, this.chunkSize, testPoint, double.PositiveInfinity, out closestPoint, out distanceSqrd);
        }

        // visit in order of bb distance
        int rootLayerItemCount = this.bbTree[0].Length;
        var order = new int[rootLayerItemCount];
        var dist = new double[rootLayerItemCount];
        for (int i = 0; i < rootLayerItemCount; i++)
        {
            order[i] = i;
            dist[i] = TreeHelper.BBDistanceToSquared(testPoint, bbTree[0][i]);
        }
        Array.Sort(dist, order);

        for (int i = 0; i < rootLayerItemCount; i++)
        {
            int idx = order[i];
            double thisDistanceSqrd;
            Point3d thisClosestPoint;
            bool inRange = this.evaluateTreeAt(0, idx, testPoint, distanceSqrd, out thisClosestPoint, out thisDistanceSqrd);
            if (!inRange) continue;
            if (thisDistanceSqrd > distanceSqrd) continue;
            distanceSqrd = thisDistanceSqrd;
            closestPoint = thisClosestPoint;
        }

        return closestPoint;
    }

    
    /**
    * Evaluate the tree at a certain depth and index for a closest point
    * this will recursively walk down this tree path until arriving at the root triangles
    * and calculate a closest point unless stopped by a negative distance inclusion check
    **/
    private bool evaluateTreeAt(
        int depth, 
        int index, 
        Point3d testPoint, 
        double maxDistanceSqrd, 
        out Point3d closestPoint, 
        out double distanceSqrd) 
    {
        if (depth < 0 || depth >= this.treeLayers) throw new Exception(OUT_OF_BOUNDS_ERROR);

        this.layerChecks[depth]++; // keep track of layer hits

        closestPoint = Point3d.Unset;
        distanceSqrd = double.PositiveInfinity;

        // check current layer for pruning
        var bbox = this.bbTree[depth][index];
        if (TreeHelper.BBDistanceToSquared(testPoint, bbox) > maxDistanceSqrd) return false;
        
        // then proceed one layer deeper
        int nextDepth = depth + 1;

        if (nextDepth < this.treeLayers) 
        {
            // collect results from deeper layer and return up the stack
            var nextLayer = this.bbTree[nextDepth];
            int offset = index * this.chunkSize;
            int itemCount = Math.Min(this.chunkSize, nextLayer.Length - offset); // guard against overshooting list end
            double best = maxDistanceSqrd;

            // visit child nodes in order of distance
            // uses small stack allocations for speed
            Span<int> order = stackalloc int[itemCount];
            Span<double> dist = stackalloc double[itemCount];
            for (int i = 0; i < itemCount; i++)
            {
                order[i] = offset + i;
                dist[i] = TreeHelper.BBDistanceToSquared(testPoint, bbTree[nextDepth][offset + i]);
            }
            dist.Sort(order);

            for (int i = 0; i < itemCount; i++) 
            {
                int nextIndex = order[i];
                double dSqrd;
                Point3d p;
                bool inRange = evaluateTreeAt(nextDepth, nextIndex, testPoint, best, out p, out dSqrd);
                if (!inRange) continue;
                if (dSqrd > best) continue;
                best = dSqrd;
                distanceSqrd = dSqrd;
                closestPoint = p;
            }
            return !double.IsInfinity(distanceSqrd);
        } 
        else 
        {
            // reached tree bottom
            // evaluate underlying triangle chunk and return closest point and distance
            int offset = index * this.chunkSize;
            bool inRange = closestPointAtTriangleChunk(
                offset, 
                this.chunkSize,
                testPoint,
                maxDistanceSqrd, 
                out closestPoint, 
                out distanceSqrd
            );
            return inRange;
        }
        
    }

    /** 
    * Iterare over triangles chunk to find closest point and distance 
    * returns true on hit, false when every triangle is outside maxDistance
    **/
    bool closestPointAtTriangleChunk(
        int startIndex, 
        int chunkSize,
        Point3d testPoint,
        double maxDistanceSqrd, 
        out Point3d closestPointOnTriangles, 
        out double closestPointDistanceSqrd
        ) 
    {
        closestPointOnTriangles = Point3d.Unset;
        closestPointDistanceSqrd = double.PositiveInfinity;

        var endIndex = Math.Min(startIndex + chunkSize, this.triangles.Length);
        for (int i = startIndex; i < endIndex; i++) 
        {
            this.triangleChecks++;

            //Triangle3d tri = this.triangles[i];
            //var pt = BVHTree.ClosestPointOnTriangle(testPoint, tri.A, tri.B, tri.C);
            var pt = TreeHelper.ClosestPointOnTriangle(testPoint, 
            this.triangleVertices[i * 3],
            this.triangleVertices[i * 3 + 1],
            this.triangleVertices[i * 3 + 2]);
            var dSqr = testPoint.DistanceToSquared(pt);
            if (dSqr > maxDistanceSqrd) continue;
            if (dSqr > closestPointDistanceSqrd) continue;
            closestPointOnTriangles = pt;
            closestPointDistanceSqrd = dSqr; 
        }

        return !double.IsInfinity(closestPointDistanceSqrd); 
    }

    static private BoundingBox BoundingBoxFromTriangles(Triangle3d[] triangles, int offset, int count)
    {
        int endIndex = Math.Min(offset + count, triangles.Length);
        Point3d[] pts = triangles
            .Skip(offset)
            .Take(endIndex - offset)
            .SelectMany(t => new[] { t.A, t.B, t.C })
            .ToArray();

        return new BoundingBox(pts);
    }

    static private BoundingBox BoundingBoxFromBoundingBoxes(BoundingBox[] boxes, int offset, int size)
    {
        var box = BoundingBox.Empty;
        int endIndex = Math.Min(offset + size, boxes.Length);
        for (int i = offset; i < endIndex; i++)
        {
            box.Union(boxes[i]);
        }
        return box;
    }

}

/** ______________________________________________________ **/
/** code below is not written by me - generated by chatGPT **/
/** ______________________________________________________ **/

static class TreeHelper {

    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    public static double BBDistanceToSquared(in Point3d p, in BoundingBox b)
    {
        var min = b.Min;
        var max = b.Max;

        double dx = 0.0;
        double px = p.X;
        if (px < min.X) dx = min.X - px;
        else if (px > max.X) dx = px - max.X;

        double dy = 0.0;
        double py = p.Y;
        if (py < min.Y) dy = min.Y - py;
        else if (py > max.Y) dy = py - max.Y;

        double dz = 0.0;
        double pz = p.Z;
        if (pz < min.Z) dz = min.Z - pz;
        else if (pz > max.Z) dz = pz - max.Z;

        return dx * dx + dy * dy + dz * dz;
    }

    
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    static double Dot(
        double ax, double ay, double az,
        double bx, double by, double bz)
        => ax * bx + ay * by + az * bz;


    /** Ericson-style closest point on triangle routine **/
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    public static Point3d ClosestPointOnTriangle(in Point3d p, in Point3d a, in Point3d b, in Point3d c)
    {
        // Load scalars once
        double px = p.X, py = p.Y, pz = p.Z;
        double ax = a.X, ay = a.Y, az = a.Z;
        double bx = b.X, by = b.Y, bz = b.Z;
        double cx = c.X, cy = c.Y, cz = c.Z;

        // ab = b - a, ac = c - a, ap = p - a
        double abx = bx - ax, aby = by - ay, abz = bz - az;
        double acx = cx - ax, acy = cy - ay, acz = cz - az;
        double apx = px - ax, apy = py - ay, apz = pz - az;

        double d1 = Dot(abx, aby, abz, apx, apy, apz);
        double d2 = Dot(acx, acy, acz, apx, apy, apz);
        if (d1 <= 0.0 && d2 <= 0.0) return a; // barycentric (1,0,0)

        // bp = p - b
        double bpx = px - bx, bpy = py - by, bpz = pz - bz;
        double d3 = Dot(abx, aby, abz, bpx, bpy, bpz);
        double d4 = Dot(acx, acy, acz, bpx, bpy, bpz);
        if (d3 >= 0.0 && d4 <= d3) return b; // (0,1,0)

        double vc = d1 * d4 - d3 * d2;
        if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
        {
            double v = d1 / (d1 - d3);
            return new Point3d(ax + v * abx, ay + v * aby, az + v * abz); // a + v*ab
        }

        // cp = p - c
        double cpx = px - cx, cpy = py - cy, cpz = pz - cz;
        double d5 = Dot(abx, aby, abz, cpx, cpy, cpz);
        double d6 = Dot(acx, acy, acz, cpx, cpy, cpz);
        if (d6 >= 0.0 && d5 <= d6) return c; // (0,0,1)

        double vb = d5 * d2 - d1 * d6;
        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
        {
            double w = d2 / (d2 - d6);
            return new Point3d(ax + w * acx, ay + w * acy, az + w * acz); // a + w*ac
        }

        double va = d3 * d6 - d5 * d4;
        if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
        {
            // edge BC
            double bcx = cx - bx, bcy = cy - by, bcz = cz - bz;
            double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return new Point3d(bx + w * bcx, by + w * bcy, bz + w * bcz);
        }

        // inside face region
        double denom = 1.0 / (va + vb + vc);
        double vFace = vb * denom;
        double wFace = vc * denom;

        return new Point3d(
            ax + abx * vFace + acx * wFace,
            ay + aby * vFace + acy * wFace,
            az + abz * vFace + acz * wFace
        );
    }

    static BoundingBox CentroidBounds(Triangle3d[] tris)
    {
        var bb = BoundingBox.Empty;
        for (int i = 0; i < tris.Length; i++)
        {
            var t = tris[i];
            var c = new Point3d(
                (t.A.X + t.B.X + t.C.X) / 3.0,
                (t.A.Y + t.B.Y + t.C.Y) / 3.0,
                (t.A.Z + t.B.Z + t.C.Z) / 3.0
            );
            bb.Union(c);
        }
        return bb;
    }

    public static void SortTrianglesByMorton(Triangle3d[] tris, int bits = 10)
    {
        var bbox = CentroidBounds(tris);
        if (!bbox.IsValid) return;

        // compute keys once (avoid recomputing during comparisons)
        var keys = new uint[tris.Length];
        for (int i = 0; i < tris.Length; i++)
        {
            var t = tris[i];
            var c = new Point3d(
                (t.A.X + t.B.X + t.C.X) / 3.0,
                (t.A.Y + t.B.Y + t.C.Y) / 3.0,
                (t.A.Z + t.B.Z + t.C.Z) / 3.0
            );
            keys[i] = MortonCode(c, bbox, bits);
        }

        // sort indices by key, then permute triangles (simple + stable enough)
        int[] idx = new int[tris.Length];
        for (int i = 0; i < idx.Length; i++) idx[i] = i;

        Array.Sort(idx, (i, j) => keys[i].CompareTo(keys[j]));

        var sorted = new Triangle3d[tris.Length];
        for (int k = 0; k < idx.Length; k++)
            sorted[k] = tris[idx[k]];

        Array.Copy(sorted, tris, tris.Length);
    }

    static uint Quantize(double v, double min, double max, int bits)
    {
        if (max <= min) return 0;

        double t = (v - min) / (max - min);
        t = Math.Clamp(t, 0.0, 1.0);

        return (uint)(t * ((1u << bits) - 1));
    }

    // Naive interleave: O(bits)
    static uint InterleaveBits(uint x, uint y, uint z, int bits)
    {
        uint code = 0;

        for (int i = 0; i < bits; i++)
        {
            uint mask = 1u << i;
            code |= ((x & mask) >> i) << (3 * i + 2);
            code |= ((y & mask) >> i) << (3 * i + 1);
            code |= ((z & mask) >> i) << (3 * i + 0);
        }

        return code;
    }

    // Simple Morton code from a point
    static uint MortonCode(
        Point3d p,
        BoundingBox bbox,
        int bits = 10)
    {
        uint x = Quantize(p.X, bbox.Min.X, bbox.Max.X, bits);
        uint y = Quantize(p.Y, bbox.Min.Y, bbox.Max.Y, bits);
        uint z = Quantize(p.Z, bbox.Min.Z, bbox.Max.Z, bits);

        return InterleaveBits(x, y, z, bits);
    }

    public static void SortPointsByMorton(
        IReadOnlyList<Point3d> pts,
        out Point3d[] sortedPts,
        out int[] sortedToOrig,
        int bits = 10)
    {
        int n = pts.Count;

        // bounds over points
        var bb = BoundingBox.Empty;
        for (int i = 0; i < n; i++) bb.Union(pts[i]);
        if (!bb.IsValid)
        {
            sortedPts = pts.ToArray();
            sortedToOrig = Enumerable.Range(0, n).ToArray();
            return;
        }

        // keys + index list
        var keys = new uint[n];
        var idx = new int[n];
        for (int i = 0; i < n; i++)
        {
            idx[i] = i;
            keys[i] = TreeHelper.MortonCode(pts[i], bb, bits);
        }

        // sort indices by key; tie-break by index for determinism
        Array.Sort(idx, (i, j) =>
        {
            int c = keys[i].CompareTo(keys[j]);
            return c != 0 ? c : i.CompareTo(j);
        });

        // build reordered points + mapping
        sortedPts = new Point3d[n];
        sortedToOrig = new int[n];
        for (int k = 0; k < n; k++)
        {
            int orig = idx[k];
            sortedPts[k] = pts[orig];
            sortedToOrig[k] = orig;
        }
    }
}
