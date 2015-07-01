namespace CGAL {

/*!
\ingroup PkgMesh2Functions

The function `lloyd_optimize_mesh_2()` is a mesh optimization process 
based on the minimization of a global energy function. 

In `lloyd_optimize_mesh_2()`, the minimized global energy may be interpreted 
as the \f$ L^1\f$-norm of the error achieved 
when the function \f$ x^2\f$ is interpolated on the mesh domain 
using a piecewise linear function which is linear 
in each cell of the Voronoi diagram of the mesh vertices. 

The optimizer `lloyd_optimize_mesh_2()` works in iterative steps. 
At each iteration, mesh vertices are moved into 
positions that bring to zero the energy gradient 
and the Delaunay triangulation is updated. 
Vertices on the mesh boundaries are not moved.

\tparam CDT is required to be or derive from `CGAL::Constrained_Delaunay_triangulation_2`,
with vertex base and face base of its underlying  `TriangulationDataStructure_2` 
respectively implementing the concepts `DelaunayMeshFaceBase_2` and `DelaunayMeshVertexBase_2`.
The argument `cdt`, passed by reference, provides the initial mesh 
and is modified by the algorithm to represent the final optimized mesh. 

\tparam PointIterator must be an iterator with value type `CGAL::Kernel::Point_2`


The function has several optional parameters which are named parameters
(we use the Boost.Parameter library). 
Therefore, when calling the function, the parameters can be provided in any order 
provided that the names of the parameters are used 
(see example at the bottom of this page). 

\cgalHeading{Named Parameters}

- <b>`parameters::time_limit`</b>
is used to set up, in seconds, 
a CPU time limit after which the optimization process is stopped. This time is 
measured using `CGAL::Timer`. 
The default value is 0 and means that there is no time limit. 
\pre `time_limit` \f$ \geq\f$ 0

- <b>`parameters::%max_iteration_number`</b> sets a limit on the 
number of performed iterations. The default value of 0 means that there is 
no limit on the number of performed iterations. 
\pre `max_iteration_number`\f$ \geq\f$ 0

- <b>`parameters::%convergence`</b> is a stopping criterion based on convergence: 
the optimization process is stopped, when at the last iteration, 
the displacement of any vertex is less than a given fraction of the 
length of the shortest edge incident to that vertex. 
The parameter `convergence` gives the threshold ratio. 
\pre 0 \f$ \leq\f$ `convergence` \f$ \leq\f$ 1

- <b>`parameters::freeze_bound`</b> is designed to reduce running time of each 
optimization iteration. Any vertex that has a displacement less than a given
fraction of the length of its shortest incident edge, is frozen (i.e.\ is 
not relocated). The parameter `freeze_bound` gives the threshold ratio.
The default value is 0.001. If it is set to 0, freezing of vertices is disabled.
\pre 0 \f$ \leq\f$ `freeze_bound` \f$ \leq\f$ 1 

- <b>`parameters::seeds_begin`</b> and <b>`parameters::seeds_end`</b>
are begin and end input iterators to iterate on seed points.
The sequence [`parameters::seeds_begin`, `parameters::seeds_end`)
defines the domain in which the mesh was generated, and should be optimized.

- <b>`parameters::mark`</b>. If `mark` is set to true, the mesh domain
is the union of the bounded connected components including at least one seed.
If `mark` is false, the domain is the union of the bounded components including
no seed. Note that the unbounded component of the plane is never optimized.
The default value is false.




\return
The function `lloyd_optimize_mesh_2()` returns a value of type `CGAL::Mesh_optimization_return_code` 
which is:
<UL> 
<LI>`CGAL::TIME_LIMIT_REACHED` when the time limit is reached. 
<LI>`CGAL::MAX_ITERATION_NUMBER_REACHED` when `lloyd_optimize_mesh_2()` stops because it has performed `max_iteration_number` iterations. 
<LI>`CGAL::CONVERGENCE_REACHED` when `lloyd_optimize_mesh_2()` stops because the convergence criterion 
is met.
<LI>`CGAL::ALL_VERTICES_FROZEN` when all vertices have been frozen, when the 
`freeze_bound` parameter is set to a positive value.
<LI>`CGAL::CANT_IMPROVE_ANYMORE` when `lloyd_optimize_mesh_2()` stops because 
most vertices have been frozen, and no better convergence can be reached.
</UL> 

\cgalHeading{Example}


\code{.cpp} 
// Lloyd-smoothing until convergence reaches 0.01, freezing vertices which 
// move less than 0.001*shortest_incident_edge_length 
lloyd_optimize_mesh_2(cdt,
                      parameters::convergence=0.01, 
                      parameters::freeze_bound=0.001); 

\endcode

\sa `Mesh_optimization_return_code` 
\sa `CGAL::refine_Delaunay_mesh_2()`

*/

template<typename CDT, typename PointIterator>
Mesh_optimization_return_code
lloyd_optimize_mesh_2(CDT& cdt,
  double parameters::time_limit=0,
  std::size_t parameters::max_iteration_number=0,
  double parameters::convergence=0.001,
  double parameters::freeze_bound = 0.001,
  PointIterator parameters::seeds_begin = PointIterator(),
  PointIterator parameters::seeds_end = PointIterator(),
  bool parameters::mark = false);

} /* namespace CGAL */
