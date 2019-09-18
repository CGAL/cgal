namespace CGAL {

/*!
\ingroup PkgPeriodic3Mesh3Functions

The function `make_periodic_3_mesh_3()` is a 3D periodic mesh generator.
It produces simplicial meshes which discretize 3D periodic domains.

The periodic mesh generation algorithm is a Delaunay refinement process
followed by an optimization phase. The criteria driving the Delaunay refinement
process may be tuned to achieve the user needs with respect to
the size of mesh elements, the accuracy of boundaries approximation,
etc.

The optimization phase is a sequence of optimization processes,
amongst the following available optimizers: an ODT-smoothing,
a Lloyd smoothing, a sliver perturber, and a sliver exuder.
Each optimization process can be activated or not, according to the user requirements
and available time.
By default, only the perturber and the exuder are activated.
Note that the benefits of the exuder will be lost if the mesh
is further refined afterward, and that ODT-smoothing, Lloyd-smoothing,
and sliver perturber should never be called after the sliver exuder.
In the case of further refinement, only the sliver exuder can be used.

The function outputs the mesh to an object which provides iterators to
traverse the resulting mesh data structure or can be written to a file
(see \ref Periodic_3_mesh_3_section_examples ).

\tparam C3T3 is required to be a model of the concept
`MeshComplex_3InTriangulation_3`. This is the return type.
The type `C3T3` is in particular required to provide a nested type
`C3T3::Triangulation` for the 3D periodic triangulation
embedding the periodic mesh. The vertex and cell base classes of the
triangulation `C3T3::Triangulation` are required to be models of the
concepts `MeshVertexBase_3` and `Periodic_3TriangulationDSVertexBase_3`, and of
the concepts `MeshCellBase_3` and `Periodic_3TriangulationDSCellBase_3`, respectively.

\tparam MD is required to be a model of the concept `Periodic_3MeshDomain_3`,
or of the refined concept `Periodic_3MeshDomainWithFeatures_3` if the domain has corners
and curve segments that need to be accurately represented in the mesh.
The argument `domain` is the sole link through which the domain
to be discretized is known by the mesh generation algorithm.

\tparam MC is required to be a model of the concept
`MeshCriteria_3`, or a model of the refined concept `MeshCriteriaWithFeatures_3` if the domain has exposed features.
The argument `criteria` of type `MC` specifies the size and shape
requirements for mesh tetrahedra and surface facets. These criteria
form the rules which drive the refinement process. All mesh elements
satisfy those criteria at the end of the refinement process.
In addition, if the domain has features, the argument
`criteria` provides a sizing field to guide the discretization
of 1-dimensional exposed features.

\cgalHeading{Named Parameters}

- <b>`features`</b> allows
the user to specify whether 0 and 1-dimensional features have to be
taken into account or not
when the domain is a model of `Periodic_3MeshDomainWithFeatures_3`.
The type `Features` of this parameter is an internal undescribed type.
The library provides functions to construct appropriate values of that type.
<UL>
<LI>\link parameters::features() `parameters::features(domain)` \endlink
sets `features` according to the domain,
i.e.\ 0 and 1-dimensional features are taken into account if `domain` is a
`Periodic_3MeshDomainWithFeatures_3`. This is the default behavior
if parameter `features` is not specified.
<LI>`parameters::no_features()` prevents the representation
of 0 and 1-dimensional features in the mesh.
</UL>

The four additional parameters are optimization parameters.
They control which optimization processes are performed
and allow the user to tune the parameters of the optimization processes.
Individual optimization parameters are not described here as they are
internal types (see instead the documentation page of each optimizer).
For each optimization algorithm, there exist two global functions
that allow to enable or disable the optimizer:

- <b>`lloyd`</b>: `parameters::lloyd()` and `parameters::no_lloyd()` are designed to
trigger or not a call to `lloyd_optimize_periodic_3_mesh_3()` function and to set the
parameters of this optimizer. If one parameter is not set, the default value of
`lloyd_optimize_periodic_3_mesh_3()` is used for this parameter.

- <b>`ODT`</b>: `parameters::odt()` and `parameters::no_odt()` are designed to
trigger or not a call to `odt_optimize_periodic_3_mesh_3` function and
to set the parameters of this optimizer.
If one parameter is not set, the default value of
`odt_optimize_periodic_3_mesh_3()` is used for this parameter.

- <b>`perturb`</b>: `parameters::perturb()` and `parameters::no_perturb()` are designed to
trigger or not a call to `perturb_periodic_3_mesh_3` function and
to set the parameters of this optimizer. If one parameter is not set, the default value of
`CGAL::perturb_periodic_3_mesh_3` is used for this parameter, except for the time bound which is set to be
equal to the refinement CPU time.

- <b>`exude`</b>: `parameters::exude()` and `parameters::no_exude()` are designed to
trigger or not a call to `exude_periodic_3_mesh_3()` function and to override to set the
parameters of this optimizer. If one parameter is not set, the default value of
`exude_periodic_3_mesh_3()` is used for this parameter, except for the time bound which is set to be
equal to the refinement CPU time.

The optimization parameters can be passed in an arbitrary order. If one parameter
is not passed, its default value is used. The default values are
`no_lloyd()`, `no_odt()`, `perturb()` and `exude()`.

Note that whatever may be the optimization processes activated,
they are always launched in the order that is a suborder
of the following (see user manual for further
 details): *ODT-smoother*, *Lloyd-smoother*, *perturb*, *exude*.

Beware that optimization of the mesh is obtained
by perturbing mesh vertices and modifying the mesh connectivity
and that this has an impact
on the strict compliance to the refinement criteria.
Though a strict compliance to mesh criteria
is guaranteed at the end of the Delaunay refinement, this may no longer be true after
some optimization processes. Also beware that the default behavior does involve some
optimization processes.

\sa `refine_periodic_3_mesh_3()`
\sa `make_mesh_3()`

\sa `parameters::features()`
\sa `parameters::no_features()`
\sa `exude_periodic_3_mesh_3()`
\sa `perturb_periodic_3_mesh_3()`
\sa `lloyd_optimize_periodic_3_mesh_3()`
\sa `odt_optimize_periodic_3_mesh_3()`
\sa `parameters::exude()`
\sa `parameters::no_exude()`
\sa `parameters::perturb()`
\sa `parameters::no_perturb()`
\sa `parameters::lloyd()`
\sa `parameters::no_lloyd()`
\sa `parameters::odt()`
\sa `parameters::no_odt()`
*/

template <class C3T3, class MD, class MC>
C3T3 make_periodic_3_mesh_3(const MD& domain,
                            const MC& criteria,
                            parameters::internal::Features_options features = parameters::features(domain),
                            parameters::internal::Lloyd_options lloyd = parameters::no_lloyd(),
                            parameters::internal::Odt_options odt = parameters::no_odt(),
                            parameters::internal::Perturb_options perturb = parameters::perturb(),
                            parameters::internal::Exude_options exude = parameters::exude());
} /* namespace CGAL */
