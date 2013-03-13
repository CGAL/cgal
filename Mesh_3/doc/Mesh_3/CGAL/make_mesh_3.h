namespace CGAL {

/*!
\ingroup PkgMesh_3Functions

The function `make_mesh_3()` is a 3D 
mesh generator. It produces simplicial meshes which discretize 
3D domains. 

The mesh generation algorithm is a Delaunay refinement process 
followed by an optimization phase. 
The criteria driving the Delaunay refinement 
process may be tuned to achieve the user needs with respect to 
the size of mesh elements, the accuracy of boundaries approximation, 
etc. 

The optimization phase is a sequence of optimization processes, 
including possibly a Lloyd smoothing, an odt-smoothing, 
a perturber and an exuder. 
Each optimization process 
can be activated or not, 
according to the user requirements 
and available time. 
By default, only the perturber and the exuder are activated. 
Note that the benefits of the exuder will be lost if the mesh 
is further refined afterward. 

The function outputs the mesh to an object which provides iterators to 
traverse the resulting mesh data structure or can be written to a file 
(see \ref Mesh_3_section_examples ). 


\tparam C3T3 is required to be a model of 
the concept 
`MeshComplex_3InTriangulation_3`. This is the return type. 
The type `C3T3` is in particular required to provide a nested type 
`C3T3::Triangulation` for the 3D triangulation 
embedding the mesh. The vertex and cell base classes of the 
triangulation `C3T3::Triangulation` are required to be models of the 
concepts `MeshVertexBase_3` and `MeshCellBase_3` 
respectively. 

\tparam MeshDomain_3 is required to be a model of 
the concept `MeshDomain_3`, or of the refined concept 
`MeshDomainWithFeatures_3` 
if the domain has corners and curve segments that need to be accurately represented in the mesh. 
The argument `domain` 
is the sole link through which the domain 
to be discretized is known by the mesh generation algorithm. 


\tparam MeshCriteria has to be a model of the concept 
`MeshCriteria_3`, or a model of the refined concept `MeshCriteriaWithFeatures_3` if the domain has exposed features. 
The argument `criteria` of 
type `MeshCriteria` specifies the 
size and shape requirements for mesh tetrahedra 
and surface facets. These criteria 
form the rules which drive the refinement process. All mesh elements 
satisfy those criteria at the end of the refinement process. 
In addition, if the domain has features, the argument 
`criteria` provides a sizing field to guide the discretization 
of 1-dimensional exposed features. 

\cgalHeading{Named Parameters}

- <b>`features`</b> allows 
the user to specify if 0 and 1-dimensional features actually have to be 
taken into account or not 
when the domain is a model of `MeshDomainWithFeatures_3`. 
The type `Features` of this parameter is an internal undescribed type. 
The library provides functions to construct appropriate values of that type. 
<UL> 
<LI>`parameters::features(domain)` sets `features` according to the domain, 
i.e.\ 0 and 1-dimensional features are taken into account if `domain` is a 
`MeshDomainWithFeatures_3`. This is the default behavior 
if parameter `features` is not specified. 
<LI>`parameters::no_features()` prevents the representation 
of 0 and 1-dimensional features in the mesh. 
</UL> 

The four additional parameters are optimization parameters. 
They control which optimization processes are performed 
and allow the user to tune the parameters of the optimization processes. 
We do not describe the types of optimization parameters as they are 
internal types. The package defines two global 
functions for each optimization parameter 
to generate appropriate value of this parameter. 

- <b>`lloyd`</b>  `parameters::lloyd()` and `parameters::no_lloyd()` are designed to 
trigger or not a call to `lloyd_optimize_mesh_3()` function and to set the 
parameters of this optimizer. If one parameter is not set, the default value of 
`lloyd_optimize_mesh_3()` is used for this parameter. 

- <b>`odt`</b> `parameters::odt()` and `parameters::no_odt()` are designed to 
trigger or not a call to `CGAL::odt_optimize_mesh_3` function and 
to set the parameters of this optimizer 
If one parameter is not set, the default value of 
`odt_optimize_mesh_3()` is used for this parameter. 

- <b>`perturb`</b> `parameters::perturb()` and `parameters::no_perturb()` are designed to 
trigger or not a call to `CGAL::perturb_mesh_3` function and 
to set the parameters of this optimizer. If one parameter is not set, the default value of 
`CGAL::perturb_mesh_3` is used for this parameter, except for the time bound which is set to be 
equal to the refinement CPU time. 

- <b>`exude`</b> `parameters::exude()` and `parameters::no_exude()` are designed to 
trigger or not a call to `exude_mesh_3()` function and to override to set the 
parameters of this optimizer. If one parameter is not set, the default value of 
`exude_mesh_3()` is used for this parameter, except for the time bound which is set to be 
equal to the refinement CPU time. 


The optimization parameters can be passed in an arbitrary order. If one parameter 
is not passed, its default value is used. The default values are 
`no_lloyd()`, `no_odt()`, `perturb()` and `exude()`. 

Note that whatever may be the optimization processes activated, 
they are always launched in the order that is a suborder 
of the following (see user manual for further details): 
*lloyd*, *odt*, *perturb*, *exude*. 

Also beware, that optimization of the mesh is obtained 
by perturbing mesh vertices and modifying the mesh connectivity 
and that this has an impact 
on the strict compliance to the refinement criteria. 
Though a strict compliance to mesh criteria 
is granted at the end of the Delaunay refinement, this may no longer be true after 
some optimization processes. Also beware that the default behavior does involve some 
optimization processes. 

\sa `refine_mesh_3()` 
\sa `parameters::features()` 
\sa `parameters::no_features()` 
\sa `exude_mesh_3()` 
\sa `perturb_mesh_3()`
\sa `lloyd_optimize_mesh_3()` 
\sa `odt_optimize_mesh_3()` 
\sa `parameters::exude()` 
\sa `parameters::no_exude()` 
\sa `parameters::perturb()` 
\sa `parameters::no_perturb()` 
\sa `parameters::lloyd()` 
\sa `parameters::no_lloyd()` 
\sa `parameters::odt()` 
\sa `parameters::no_odt()` 

*/

template <class C3T3,
class MeshDomain_3,
class MeshCriteria>
C3T3 make_mesh_3(MeshDomain_3 domain,
MeshCriteria criteria,
parameters::internal::Features_options features = parameters::features(domain),
parameters::internal::Lloyd_options lloyd = parameters::no_lloyd(),
parameters::internal::Odt_options odt = parameters::no_odt(),
parameters::internal::Perturb_options perturb = parameters::perturb(),
parameters::internal::Exude_options exude = parameters::exude()); 

} /* namespace CGAL */
