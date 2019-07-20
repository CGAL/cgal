namespace CGAL {

/*!
\ingroup PkgPeriodic3Mesh3Functions

The function `lloyd_optimize_periodic_3_mesh_3()` is a periodic mesh optimization
process based on the minimization of a global energy function.

This function directly calls `lloyd_optimize_mesh_3()`, but is provided for convenience.
Further information can be found on the documentation of the function `lloyd_optimize_mesh_3()`.

\note This function requires the \ref thirdpartyEigen library.
*/
template<typename C3T3, typename MD>
CGAL::Mesh_optimization_return_code
lloyd_optimize_periodic_3_mesh_3(C3T3& c3t3,
                                 const MD& domain,
                                 double parameters::time_limit=0,
                                 std::size_t parameters::max_iteration_number=0,
                                 double parameters::convergence=0.02,
                                 double parameters::freeze_bound = 0.01,
                                 bool parameters::do_freeze=true);


/*!
\ingroup PkgPeriodic3Mesh3Functions

The function `odt_optimize_periodic_3_mesh_3()` is a periodic mesh optimization
process based on the minimization of a global energy function.

This function directly calls `odt_optimize_mesh_3()`, but is provided for convenience.
Further information can be found on the documentation of the function `odt_optimize_mesh_3()`.
*/
template<typename C3T3, typename MD>
Mesh_optimization_return_code
odt_optimize_periodic_3_mesh_3(C3T3& c3t3,
                               const MD& domain,
                               double parameters::time_limit=0,
                               std::size_t parameters::max_iteration_number=0,
                               double parameters::convergence=0.02,
                               double parameters::freeze_bound = 0.01,
                               bool parameters::do_freeze=true);

/*!
\ingroup PkgPeriodic3Mesh3Functions

The function `perturb_periodic_3_mesh_3()` is a mesh optimizer that
improves the quality of a Delaunay mesh by changing the positions of some vertices of the mesh.

This function directly calls `perturb_mesh_3()`, but is provided for convenience.
Further information can be found on the documentation of the function `perturb_mesh_3()`.
*/
template<typename C3T3, typename MD>
CGAL::Mesh_optimization_return_code
perturb_periodic_3_mesh_3(C3T3& c3t3,
                          const MD& domain,
                          double parameters::time_limit=0,
                          double parameters::sliver_bound=0);

/*!
\ingroup PkgPeriodic3Mesh3Functions

The function `exude_periodic_3_mesh_3()` performs a sliver exudation process
on a periodic Delaunay mesh.

The sliver exudation process consists in optimizing the weights of vertices
of the periodic weighted Delaunay triangulation in such a way that slivers disappear
and the quality of the mesh improves.

\warning This optimizer modifies the weight of vertices of the periodic triangulation
and, if called, must be the last optimizer to be called. If the mesh is refined after
this optimization has been performed, all improvements will be lost.

This function directly calls `exude_mesh_3()`, but is provided for convenience.
Further information can be found on the documentation of the function `exude_mesh_3()`.
*/
template<typename C3T3>
CGAL::Mesh_optimization_return_code
exude_periodic_3_mesh_3(C3T3& c3t3,
                        double parameters::time_limit=0,
                        double parameters::sliver_bound=0);

} /* namespace CGAL */
