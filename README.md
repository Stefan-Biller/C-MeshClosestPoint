# Mesh Closest Point (Rhino/Grasshopper)

Compute closest points from a set of input points to a mesh efficiently using an R-Tree of axis aligned bounding boxes with Morton ordering and parallel processing.

## Overview

- BVH over triangle chunks: hierarchical `BoundingBox` layers from coarse to fine.
- Morton (Z-order) sorting: reorders triangles and input points to improve locality and cache coherence.
- Parallel batches: evaluates point queries in parallel for throughput.
- Early pruning: limits search radius using the last closest point to avoid unnecessary tree traversals.
- Instrumentation: reports setup and execution timing per iteration and totals.

See the implementation in [C# MeshClosestPoint.cs](C%23%20MeshClosestPoint.cs).

## Grasshopper Usage

1. Add a C# Script component to your Grasshopper canvas.
2. Open the component, replace its template code with the contents of [C# MeshClosestPoint.cs](C%23%20MeshClosestPoint.cs).
3. Create inputs and outputs as listed below (exact names and types):

### Inputs

- `M` (Mesh): Target mesh. Quads are converted to triangles internally.
- `Pts` (List<Point3d>): Points to project onto the mesh.
- `layers` (int): Number of BVH layers (≥ 1). Typical: 3–6.
- `nodeSize` (int): Items per node/chunk (≥ 1). Typical: 16–64.
- `iterations` (int): Repeat count to profile performance. Set to 1 for normal use.

### Outputs

- `P` (IEnumerable<GH_Point>): Closest points on the mesh for the input `Pts`.
- `T` (BoundingBox[][]): BVH tree as layered bounding boxes (coarsest → finest).
- `t` (double): Total time in milliseconds of the last iteration.

The component also prints timing breakdowns and diagnostic counters to the Grasshopper console.

## Parameters and Tuning

- **`layers`**: More layers increase pruning potential but add construction overhead. 3–6 usually balances well.
- **`nodeSize`**: Larger chunks reduce tree depth and bbox count; smaller chunks increase precision near leaves. Start at 32.
- **Parallel batch size**: Internal constant `PARALLEL_BATCH_SIZE = 32`. Adjust if you change hardware or point counts substantially.
- **Tolerance**: BVH uses a small tolerance (`1e-9`) consistent with Rhino millimeter scale.

## How It Works

- Triangulation: Mesh quads are converted to triangles; a flat array of triangle vertices is built for tight loops.
- BVH Build: Leaf layer contains bounding boxes over triangle chunks; upper layers union child boxes.
- Morton Ordering: Triangles are sorted by centroid Morton codes; points are also sorted by Morton to improve traversal locality.
- Closest Point Query:
  - Computes initial distance from a small triangle sample if no bound exists.
  - Traverses BVH layers in distance order, pruning subtrees whose bbox distance exceeds the current best.
  - Evaluates triangle chunks using an Ericson-style closest point routine.
- Parallel Execution: Points are processed in fixed-size batches using `Parallel.ForEach`.

## Console Diagnostics

The script prints per-iteration and aggregated totals:

- Setup: BVH build and point sorting time.
- Execution: Query time for all points.
- Tree layer checks: Visits per layer for the last run.
- Triangle checks: Number of triangle evaluations; per-point average.
- Average time across iterations.

Use these to tune `layers` and `nodeSize` for your meshes and point distributions.

## Limitations

- Requires a valid mesh convertible to triangles; throws if triangulation fails.
- `T` (BVH) is a jagged array of `BoundingBox`; visualization requires custom drawing if you want to see the tree.
- Memory footprint depends on mesh size and `nodeSize`.

## Tips

- If points are spatially coherent, Morton sorting and the adaptive radius greatly reduce work.
- For highly scattered points, favor slightly larger `nodeSize` to reduce bbox overhead.
- For very dense meshes, increase `layers` to improve pruning depth.

## Attribution

Portions of helper routines in `TreeHelper` (Morton code and triangle closest point) are noted as generated in-code. The overall design and integration are tailored for Grasshopper + RhinoCommon.

## Path and Files

- Main script: [C# MeshClosestPoint.cs](C%23%20MeshClosestPoint.cs)
- README: [README.md](README.md)
