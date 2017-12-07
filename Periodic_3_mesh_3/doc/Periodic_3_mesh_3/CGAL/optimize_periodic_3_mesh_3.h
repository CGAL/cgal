namespace CGAL {

/*!
\ingroup PkgPeriodic_3_mesh_3Functions

The function `lloyd_optimize_periodic_3_mesh_3()` is a periodic mesh optimization
process based on the minimization of a global energy function.

\pre `time_limit` \f$ \geq\f$ 0 and 0 \f$ \leq\f$ `convergence` \f$ \leq\f$ 1 and 0 \f$ \leq\f$ `freeze_bound` \f$ \leq\f$ 1

\tparam MD is required to be a model of the concept `Periodic_3MeshDomain_3`,
or of the refined concept `Periodic_3MeshDomainWithFeatures_3` if the domain has corners
and curve segments that need to be accurately represented in the mesh.

This function directly calls `lloyd_optimize_mesh_3()`, but is provided for convience.

*/
template<typename C3T3, typename MD>
CGAL::Mesh_optimization_return_code
lloyd_optimize_periodic_3_mesh_3(C3T3& c3t3,
                                 MD domain,
                                 double parameters::time_limit=0,
                                 std::size_t parameters::max_iteration_number=0,
                                 double parameters::convergence=0.02,
                                 double parameters::freeze_bound = 0.01,
                                 bool parameters::do_freeze=true);


/*!
\ingroup PkgPeriodic_3_mesh_3Functions

The function `odt_optimize_periodic_3_mesh_3()` is a periodic mesh optimization
process based on the minimization of a global energy function.

\pre `time_limit` \f$ \geq\f$ 0 and 0 \f$ \leq\f$ `convergence` \f$ \leq\f$ 1 and 0 \f$ \leq\f$ `freeze_bound` \f$ \leq\f$ 1

\tparam MD is required to be a model of the concept `Periodic_3MeshDomain_3`,
or of the refined concept `Periodic_3MeshDomainWithFeatures_3` if the domain has corners
and curve segments that need to be accurately represented in the mesh.

This function directly calls `odt_optimize_mesh_3()`, but is provided for convience.

*/
template<typename C3T3, typename MD>
Mesh_optimization_return_code
odt_optimize_periodic_3_mesh_3(C3T3& c3t3,
                               MD domain,
                               double parameters::time_limit=0,
                               std::size_t parameters::max_iteration_number=0,
                               double parameters::convergence=0.02,
                               double parameters::freeze_bound = 0.01,
                               bool parameters::do_freeze=true);

/*!
\ingroup PkgPeriodic_3_mesh_3Functions

The function `perturb_periodic_3_mesh_3()` is a mesh optimizer that
improves the quality of a Delaunay mesh by changing the mesh vertices positions.

\pre `time_limit` \f$ \geq\f$ 0 and 0 \f$ \leq\f$ `sliver_bound` \f$ \leq\f$ 180

\tparam MD is required to be a model of the concept `Periodic_3MeshDomain_3`,
or of the refined concept `Periodic_3MeshDomainWithFeatures_3` if the domain has corners
and curve segments that need to be accurately represented in the mesh.

This function directly calls `perturb_mesh_3()`, but is provided for convience.

*/
template<typename C3T3, typename MD>
CGAL::Mesh_optimization_return_code
perturb_periodic_3_mesh_3(C3T3& c3t3,
                          MD domain,
                          double parameters::time_limit=0,
                          double parameters::sliver_bound=0);

/*!
\ingroup PkgPeriodic_3_mesh_3Functions

The function `exude_periodic_3_mesh_3()` performs a sliver exudation process
on a periodic Delaunay mesh.

The sliver exudation process consists in turning the periodic Delaunay triangulation
into a periodic weighted Delaunay triangulation and optimizing the weights
of vertices in such a way that slivers disappear and the quality of the mesh improves.

\pre `time_limit` \f$ \geq\f$ 0 and 0 \f$ \leq\f$ `sliver_bound` \f$ \leq\f$ 180

\tparam MD is required to be a model of the concept `Periodic_3MeshDomain_3`,
or of the refined concept `Periodic_3MeshDomainWithFeatures_3` if the domain has corners
and curve segments that need to be accurately represented in the mesh.

This function directly calls `exude_mesh_3()`, but is provided for convience.

*/
template<typename C3T3>
CGAL::Mesh_optimization_return_code
exude_periodic_3_mesh_3(C3T3& c3t3,
                        double parameters::time_limit=0,
                        double parameters::sliver_bound=0);

} /* namespace CGAL */
