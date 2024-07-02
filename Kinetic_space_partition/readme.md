** GENERAL **

All the test data are placed in the folder examples/data. The name of the folder indicates the type of data or its complexity.
The complexity 1 is low while the complexity 6 is high.

Examples are the entry point to the package:
- kinetic_2d_example: 2D kinetic algorithm on a random set of segments
- kinetic_precomputed_shapes_example: 3D kinetic algorithm on a set of user-defined polygons read from a file
- kinetic_random_shapes_example: 3D kinetic algorithm on a set of random polygons, you can provide the number
  of input polygons to be generated and the number of vertices each polygon should have
- kinetic_reconstruction_example: given a ply with classified point cloud, it first fits polygons using Region Growing,
  then runs 3D kinetic, and finally filters out exterior volumes creating a reconstructed model

Tests:
- kinetic_2d_stress_test: tests the 2D kinetic algorithm
- kinetic_3d_test_all: tests the 3D kinetic algorithm on all data from the examples folder

The 3D algorithm is running on macOS with zero warnings. On Windows and Linux, it compiles and works,
but can generate warnings related to conversions between std:size_t and int e.g. or similar warnings.
The 2D algorithm should work on all platforms as well.


** CURRENT ISSUES **


-- EPICK versus EPECK --

By running:
./kinetic_precomputed_shapes_example data/stress-test-4/test-9-rnd-polygons-12-4.off

first with EPICK and then with EPECK shows a huge difference in runtime. This
data set contains 12 random polygons, each having 4 vertices. The time for EPICK is
milliseconds while for EPECK is about 3 minutes. It can also be noticed that most of
the time is spent at the last iterations while the first iterations are very fast.
It is even slower for `Simple_cartesian<Gmpq>`.

It is probably happening because at the latest iterations the computation involves
all operations carried out from the first iteration. Evaluating these cascading operations
is a slow process. E.g. computing an intersection between two lines at the first iterations
takes seconds while at the higher iterations it takes minutes.

Another issue shows up when running the code with the following data set:
./kinetic_precomputed_shapes_example data/real-data-test/test-40-polygons.ply 6
that is with k = 6, where k is the number of allowed intersections between support planes
of the input polygons. Here, at about 5300 iteration, assertions start failing but not because
some parts of the algorithm are not implemented (as indicated in the TODO) but because
the scaling is no longer uniform due to the accumulating errors. It leads to wrong configurations
and hence to the wrong events, which in practice with the correct precision cannot happen.


-- Possible solutions --

- Use EPICK and compute always with respect to the input. E.g. intersections between lines should be
done not between lines at the current and previous iterations but between lines at the current
and first iterations. The reason for this optimization is that if we use simply EPICK,
when the number of input polygons grow, we bump into an issue of the accumulating error that gets
bigger and bigger with each iteration and at some point can break the results. It does not happen
with EPECK but we lose speed there.

- Use EPECK only when it is absolutely necessary like when computing intersections and
use EPICK for all other computations. This way we avoid accumulating errors and should
preserve speed.

A few ideas shortly after the meeting with Sebastien:
- speed up kinetic with EPECK by optimizing the parts, which are slow, vtune it
- use different speed for each polygon edge during the uniform scaling, a similar
  trick to what is used in the straight skeleton now
- point_2() wrapper inside Support_plane.h: it can compute point not based on the point constructed
  at the previous event (which is cascading) but based on the interesting polygon edge and intersection graph edge
- avoid certain events in the queue or keep track only of the last 10 events instead of all of them


** INTERNALS **


-- File descriptions --

- KSR sub-folder, used in both KSR2 and KSR3:

  - conversions.h: the kinetic_traits class that performs intersections, all intersections
  in the code are called from this class, here the hybrid mode could be potentially implemented.

  - debug.h: a file with a bunch of functions which enable to export / dump to a file
  different intermediate steps, events, data structure, etc.

  - enum.h: enumerations used by the reconstruction code from the Reconstruction.h.

  - parameters.h: internal parameters used in the code.

  - property_map.h: property maps used by the reconstruction code from the Reconstruction.h.

  - utils.h: different internal utilities such as normals, distances, angles, etc.
  the most important one is tolerance(): this is a global tolerance used everywhere in the code

- KSR_2 sub-folder, 2D kinetic algorithm, works and tested, fully implemented by Simon.

- KSR_3 sub-folder, 3D kinetic algorithm + 3D reconstruction algorithm:

  - Data_structure.h: a class with all conversions and all basic operations, which are
  performed on the polygons and during events, it also stores the main intersection graph
  and all support planes.

  - Event_queue.h: a wrapper around boost multi_index_container that represents a queue of
  all events, which can be sorted by time or vertex type. It is based on doubles so even
  when we use EPECK globally, we need to convert to double here because boost multi_index_container
  fails to work with EPECK due to the fact that values represented by EPECK can vary when converted
  to_double() or to_interval().

  - Event.h: represents different event types used in the propagation.

  - Intersection_graph.h: a boost adjacency_list graph that stores an initial intersection graph
  fully in 3D that is used as constraints for kinetic propagation.

  - Support_plane.h: a wrapper around a support plane for each input polygon that stores a 2D
  surface mesh that represents an input polygon. So, initially this mesh has only 1 face that is
  input polygon while during propagation it is updated and new vertices, faces are added. We preserve
  the uniform scaling of the input polygon and validity of this mesh after each event. If any of these
  two fails, the whole algorithm fails. At the end of the support plane header, there is an overload
  of the operator()== that is used to compare if two planes are equal. It heavily depends on the
  tolerance parameters and in case two planes are found to be equal, they are merged at the initialization step.

  - Initializer.h: a class that gets input polygons, creates an enlarged bounding box, removes all
  near-collinear vertices of each polygon, merges all near-coplanar support planes, then inserts each
  support plane and intersects it with the bounding box. At the final step, it calls Polygon_splitter
  class in order to intersect all polygons within the bounding box, one support plane at a time.

  - Polygon_splitter.h: given a support plane, its original polygon, and a set of intersection edges created
  by intersecting all support planes with the bounding box (see Initializer.h), it builds a CDT, inserts all these edges
  in the CDT, marks all interior and exterior faces and creates a proper 2D surface mesh. This class can be parameterized
  by any kernel independently of the input kernel.

  - Propagation.h: is called after Initializer.h has done its work, this class creates an initial queue of
  events for all vertices of all polygons and handles these events one by one. This is the most time-consuming
  part of the total algorithm. The queue is updated and rebuilt until no valid events show up. This class can be parameterized
  by any kernel independently of the input kernel.

  - Finalizer.h: is called after Propagation.h, this class receives a set of 2D surface meshes from the propagation
  and an initial intersection graph and does three things: first it checks the validity of the results, it then searches
  for dangling / hanging faces if any and fills the gaps; and finally creates 3D volumes bounded by all faces of the 2D surface meshes.
  This class can be parameterized by any kernel independently of the input kernel.

  - Reconstruction.h: this class first gets a point cloud with or without semantic labels. If the labels are missing, some
  parts of the code have to be still finished. However, if labels are present, it assumes they come from urban areas that is they are
  walls, roofs, ground, and trees. It removes all tree points, fits a polygon to the ground points, runs region growing
  separately for wall and roof points and detects convex planar polygons, which approximate the input point cloud.
  Next, it runs 3D kinetic on these input polygons + ground polygon and gets a set of volumes, which partition the 3D space.
  It then estimates a probability of each volume to be inside or outside the point cloud using Visibility.h and runs the graph cut
  algorithm on these estimations using Graphcut.h. At the end, it removes all volumes, which have been labeled as exterior volumes
  and returns a reconstructed model: 3D surface mesh that approximates the input point cloud (or buildings in case labels are urban related).

  - Visibility.h: estimates a probability of each partition volume to be inside or outside a point cloud by using normals
  of the input points and sampling each volume.

  - Graphcut.h: takes initially estimated probabilities for each volume to be interior or exterior with respect to the point cloud,
  computed by Visibility.h, creates a graph where each node is a volume and each edge connects volumes to their neighboring volumes,
  and runs the mincut - maxflow Boykov optimization algorithm to define which volumes are inside and outside the point cloud.
  All edges are weighted by the face area that adjacent to two incident volumes and all nodes are weighted by the volume itself.

- Kinetic_shape_reconstruction_2.h: is an entry point to the 2D kinetic propagation.

- Kinetic_shape_reconstruction_3.h: is an entry point to the 3D kinetic propagation and 3D kinetic reconstruction algorithms.


-- Epsilon usage inside the code --

Note that epsilon tolerance is used throughout the code. It is defined inside utils.h
in the function tolerance(). It is used everywhere where we expect a result of certain precision
but it may not be the case. We check if we outside this tolerance and apply different sub-algorithms
in order to be sure we that will generate the right results. This is mostly used for EPICK. When using EPECK,
the precision is higher and the tolerance should be satisfied until there is a bug.


-- Important parts of the code --

  - point_2() wrapper from the Support_plane.h. Currently all original points are 3D,
  this wrapper takes a point, and constructs the corresponding 2D point for the chosen
  2D support plane and at the chosen time step. That means all these points are constructed
  that is a weak point when using EPECK because it is slow.

  - operator==() at the bottom of the Support_plane.h: controls if two near-coplanar planes should
  be merged. If we merge too many planes because our input parameters are very weak, we fail to create
  a valid partition. If we do not merge at all, the near-coplanar support planes may lead to intricated
  errors in the kinetic propagation due to numerical instabilities.

  - parameters: all all explained in the parameters.h. The parameters used in the Reconstruction.h are defined
  in the file examples/include/Parameters.h.

  - FUTURE POINTS AND DIRECTIONS section at the bottom of the Data_structure.h, this is where the future points
  are computed and this is a part of the code that leads to multiple precision issues, identifying a future
  point from the previous event is hard, so instead we simply translate the lines and intersect them at the next
  time step to get the point at which our current vertex will be later, but these intersections are imprecise.
  If we lose precision here, we fail to scale polygon uniformly so our point can end up behind the bounding box
  or in the direction opposite to the one we need, especially if the lines that we intersect are near parallel.
