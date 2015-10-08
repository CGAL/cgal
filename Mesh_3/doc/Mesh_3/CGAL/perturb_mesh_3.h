namespace CGAL {

/*!
\ingroup PkgMesh_3Functions

The function `perturb_mesh_3()` is a mesh optimizer that 
improves the quality of a Delaunay mesh 
by changing the mesh vertices positions. 

The perturber tries to improve the dihedral angles of the worst cells in the mesh 
degree by degree: the 
step number `n` is considered as successful 
if after this step the worst tetrahedron of the mesh has a minimal dihedral 
angle larger than `n` degrees. 
The perturber exits if this is not the case. 

\pre `time_limit` \f$ \geq\f$ 0 and 0 \f$ \leq\f$ `sliver_bound` \f$ \leq\f$ 180 

\tparam C3T3 is required to be a model of the concept 
`MeshComplex_3InTriangulation_3`. 
The argument `c3t3`, passed by 
reference, provides the initial mesh 
and is modified by the algorithm 
to represent the final optimized mesh. 

\tparam MeshDomain_3 is required to be a model of the concept 
`MeshDomain_3`. The argument `domain` must be the `MeshDomain_3` 
object used to create the `c3t3` parameter. 

The function has two optional parameters which are named parameters (we use the Boost.Parameter library). 
Therefore, when calling the function, the parameters can be provided in any order 
provided that the names of the parameters are used 
(see example at the bottom of this page). 

\cgalHeading{Named Parameters}

- <b>`parameters::time_limit`</b>
is used to set up, in seconds, 
a CPU time limit after which the optimization process is stopped. This time is 
measured using `Real_timer`. 
The default value is 0 and means that there is no time limit. 

- <b>`parameters::sliver_bound`</b>
is designed to give, in degree, a targeted 
lower bound on dihedral angles of mesh cells. 
The function `perturb_mesh_3()` runs as long as steps are successful 
and step number `sliver_bound` (after which 
the worst tetrahedron in the mesh has a smallest angle larger than 
`sliver_bound` degrees) has not been reached. 
The default value is 0 and means that there is no targeted bound: 
the perturber then runs as long as 
steps are successful. 


\return
The function `perturb_mesh_3()` returns a value of type `CGAL::Mesh_optimization_return_code` 
which is: 
<UL> 
<LI>`CGAL::BOUND_REACHED` when the targeted bound for the smallest dihedral angle in the mesh is reached. 
<LI>`CGAL::TIME_LIMIT_REACHED` when the time limit is reached. 
<LI>`CGAL::CANT_IMPROVE_ANYMORE` when the perturbation process stops because the last step is unsuccessful. 
</UL> 


\cgalHeading{Example}

\code{.cpp} 
// Perturb until every dihedral angle of the mesh is >= 10 degrees 
// No time bound is set 
perturb_mesh_3(c3t3, 
               domain, 
               parameters::sliver_bound = 10); 
\endcode 

\sa `CGAL::Mesh_optimization_return_code` 
\sa `CGAL::make_mesh_3()` 
\sa `CGAL::refine_mesh_3()` 
\sa `CGAL::exude_mesh_3()` 
\sa `CGAL::lloyd_optimize_mesh_3()` 
\sa `CGAL::odt_optimize_mesh_3()` 

*/

template<typename C3T3, typename MeshDomain_3>
Mesh_optimization_return_code
perturb_mesh_3(C3T3& c3t3,
MeshDomain_3 domain,
double parameters::time_limit=0,
double parameters::sliver_bound=0);

} /* namespace CGAL */
