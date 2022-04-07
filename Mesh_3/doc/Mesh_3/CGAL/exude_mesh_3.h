namespace CGAL {

/*!
\ingroup PkgMesh3Functions

The function `exude_mesh_3()` performs a sliver exudation process on a Delaunay mesh.

The sliver exudation process consists in optimizing the weights of vertices
of the weighted Delaunay triangulation in such a way that slivers disappear and
the quality of the mesh improves.

\warning This optimizer modifies the weight of vertices of the triangulation and,
if called, must be the last optimizer to be called. If the mesh is refined after
this optimization has been performed, all improvements will be lost.

\pre `time_limit` \f$ \geq\f$ 0 and 0 \f$ \leq\f$ `sliver_bound` \f$ \leq\f$ 180

\tparam  C3T3 is required to be a model of the concept
`MeshComplex_3InTriangulation_3`.
The argument `c3t3`, passed by
reference, provides the initial mesh
and is modified by the algorithm
to represent the final optimized mesh.

The function has two optional parameters which are named parameters (we use the Boost.Parameter library).
Therefore, when calling the function, the parameters can be provided in any order
provided that the names of the parameters are used
(see example at the bottom of this page).

\cgalHeading{Named Parameters}
- <b>`parameters::time_limit`</b> is used to set up, in seconds,
a CPU time limit after which the optimization process is stopped. This time is
measured using the `Real_timer` class.
The default value is 0 and means that there is no time limit.

- <b>`parameters::sliver_bound`</b> is designed to give, in degrees, a targeted
lower bound on dihedral angles of mesh cells.
The exudation process considers in turn all the mesh cells
that have a smallest dihedral angle less than `sliver_bound`
and tries to make them disappear by weighting their vertices.
The optimization process
stops when every cell in the mesh achieves this quality.
The default value is 0 and means that there is no targeted bound:
the exuder then runs as long as
it can improve the smallest dihedral angles of the set of cells
incident to some vertices.

\return
The function `exude_mesh_3()` returns a value of type `CGAL::Mesh_optimization_return_code`
which is:
<UL>
<LI>`CGAL::BOUND_REACHED` when the targeted bound for the smallest dihedral angle in the mesh is reached.
<LI>`CGAL::TIME_LIMIT_REACHED` when the time limit is reached.
<LI>`CGAL::CANT_IMPROVE_ANYMORE` when exudation process stops because it can no longer improve
the smallest dihedral angle of the set of cells incident to some vertex in the mesh.
</UL>

\cgalHeading{Example}

\code{.cpp}
// Exude without sliver_bound, using at most 10s CPU time
exude_mesh_3(c3t3,
             parameters::time_limit=10);
\endcode

\sa `CGAL::Mesh_optimization_return_code`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`
\sa `CGAL::perturb_mesh_3()`
\sa `CGAL::lloyd_optimize_mesh_3()`
\sa `CGAL::odt_optimize_mesh_3()`

*/

template<typename C3T3>
Mesh_optimization_return_code
exude_mesh_3(C3T3& c3t3,
             double parameters::time_limit=0,
             double parameters::sliver_bound=0);

} /* namespace CGAL */
