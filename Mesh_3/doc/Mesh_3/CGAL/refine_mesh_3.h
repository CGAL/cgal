namespace CGAL {

/*!
\ingroup PkgMesh3Functions

The function `refine_mesh_3()` is a 3D
mesh generator. It produces simplicial meshes which discretize
3D domains.

The mesh generation algorithm is a Delaunay refinement process
followed by an optimization phase.
The criteria driving the Delaunay refinement
process may be tuned to achieve the user needs with respect to
the size of mesh elements, the accuracy of boundaries approximation,
etc.

The optimization phase is a sequence of optimization processes,
amongst the following available optimizers: an ODT-smoothing,
a Lloyd smoothing, a sliver perturber, and a sliver exuder.
Each optimization process
can be activated or not,
according to the user requirements
and available time.
By default, only the perturber and the exuder are activated.
Note that the benefits of the exuder will be lost if the mesh
is further refined afterward.

\attention The function template `refine_mesh_3()` may be used to refine a previously
computed mesh, e.g.:
\code{.cpp}
C3T3 c3t3 = CGAL::make_mesh_3<C3T3>(domain,criteria);

CGAL::refine_mesh_3(c3t3, domain, new_criteria);
\endcode

Please note that we guarantee the result if and only if the domain does
not change from one refinement to the next one.


\tparam  C3T3 is required to be a model of
the concept
`MeshComplex_3InTriangulation_3`.
The argument `c3t3` is passed by
reference as this object is modified by the refinement process. As the
refinement process only adds points to the triangulation, all
vertices of the triangulation of `c3t3` remain in the
mesh during the refinement process. Object `c3t3` can be used to insert
specific points in the domain to ensure that they will be contained in the
final triangulation.
The type `C3T3` is in particular required to provide a nested type
`C3T3::Triangulation` for the 3D triangulation
embedding the mesh. The vertex and cell base classes of the
triangulation `C3T3::Triangulation` are required to be models of the
concepts `MeshVertexBase_3` and `MeshCellBase_3`
respectively.

\tparam MD is required to be a model of
the concept `MeshDomain_3` or of the refined concept
`MeshDomainWithFeatures_3` if 0 and 1-dimensional features
of the input complex have to be accurately represented in the mesh.
The argument `domain`
is the sole link through which the domain
to be discretized is known by the mesh generation algorithm.

\tparam MC is required to be a model of the concept
`MeshCriteria_3`, or a model of the refined concept `MeshCriteriaWithFeatures_3`
if the domain has exposed features. The argument `criteria` of
type `MC` specifies the
size and shape requirements for mesh tetrahedra
and surface facets. These criteria
form the rules which drive the refinement process. All mesh elements
satisfy those criteria at the end of the refinement process.
In addition, if the domain has features, the argument
`criteria` provides a sizing field to guide the discretization
of 1-dimensional exposed features.

The four additional parameters are optimization parameters.
They control which optimization processes are performed
and allow the user to tune the parameters of the optimization processes.
Individual optimization parameters are not described here as they are
internal types (see instead the documentation page of each optimizer).
For each optimization algorithm, there exist two global functions
that allow to enable or disable the optimizer:

\cgalHeading{Named Parameters}
- <b>`manifold`</b> allows the user to drive the meshing algorithm,
and ensure that the output mesh surface follows the given manifold criterion.
It can be activated with `parameters::manifold()`, `parameters::manifold_with_boundary()`
and `parameters::non_manifold()`. Note that the meshing algorithm cannot generate a manifold
surface if the input surface is not manifold.

- <b>`lloyd`</b>  `parameters::lloyd()` and `parameters::no_lloyd()` are designed to
trigger or not a call to `lloyd_optimize_mesh_3()` function and to set the
parameters of this optimizer. If one parameter is not set, the default value of
`lloyd_optimize_mesh_3()` is used for this parameter.

- <b>`ODT`</b> `parameters::odt()` and `parameters::no_odt()` are designed to
trigger or not a call to `odt_optimize_mesh_3()` function and
to set the parameters of this optimizer.
If one parameter is not set, the default value of
`odt_optimize_mesh_3()` is used for this parameter.

- <b>`perturb`</b> `parameters::perturb()` and `parameters::no_perturb()` are designed to
trigger or not a call to `perturb_mesh_3()` function and
to set the parameters of this optimizer. If one parameter is not set, the default value of
`perturb_mesh_3()` is used for this parameter, except for the time bound which is set to be
equal to the refinement CPU time.

- <b>`exude`</b> `parameters::exude()` and `parameters::no_exude()` are designed to
trigger or not a call to `exude_mesh_3()` function and to override to set the
parameters of this optimizer. If one parameter is not set, the default value of
`exude_mesh_3()` is used for this parameter, except for the time bound which is set to be
equal to the refinement CPU time.

The optimization parameters can be passed in arbitrary order. If one parameter
is not passed, its default value is used. The default values are
`no_lloyd()`, `no_odt()`, `perturb()` and `exude()`.
Note that whatever may be the optimization processes activated,
they are always launched in the order that is a suborder
of the following (see user manual for further
details): *ODT-smoother*, *Lloyd-smoother*, *perturber*, and *exuder*.

Beware that optimization of the mesh is obtained
by perturbing mesh vertices and modifying the mesh connectivity
and that this has an impact
on the strict compliance to the refinement criteria.
Though a strict compliance to mesh criteria
is guaranteed at the end of the Delaunay refinement, this may no longer be true after
some optimization processes. Also beware that the default behavior does involve some
optimization processes.

\sa `CGAL::make_mesh_3()`
\sa `CGAL::parameters::manifold`
\sa `CGAL::parameters::manifold_with_boundary`
\sa `CGAL::parameters::non_manifold`
\sa `CGAL::exude_mesh_3()`
\sa `CGAL::perturb_mesh_3()`
\sa `CGAL::lloyd_optimize_mesh_3()`
\sa `CGAL::odt_optimize_mesh_3()`
\sa `CGAL::parameters::exude`
\sa `CGAL::parameters::no_exude`
\sa `CGAL::parameters::perturb`
\sa `CGAL::parameters::no_perturb`
\sa `CGAL::parameters::lloyd`
\sa `CGAL::parameters::no_lloyd`
\sa `CGAL::parameters::odt`
\sa `CGAL::parameters::no_odt`

*/

template <class C3T3, class MD, class MC>
void refine_mesh_3(C3T3& c3t3,
                   const MD& mesh_domain,
                   const MC& mesh_criteria,
                   parameters::internal::Lloyd_options lloyd = parameters::no_lloyd(),
                   parameters::internal::Odt_options odt = parameters::no_odt(),
                   parameters::internal::Perturb_options perturb = parameters::perturb(),
                   parameters::internal::Exude_options exude = parameters::exude(),
                   parameters::internal::Manifold_options manifold = parameters::non_manifold());

namespace parameters {

  /*!
  \ingroup PkgMesh3Parameters

  The function `parameters::manifold()` is used to drive the
  meshing algorithm for surfaces.
  It ensures that the surface of the output mesh is a manifold surface
    without boundaries.
  The manifold property of the output mesh can be achieved only if the input surface
  is a manifold.
  Note that the meshing algorithm provably terminates only if the input
  sharp edges have been protected, using the
  feature protection (see \ref Mesh_3Protectionof0and1dimensionalExposed).

  \sa `CGAL::make_mesh_3()`
  \sa `CGAL::refine_mesh_3()`
  \sa `CGAL::parameters::manifold_with_boundary()`
  \sa `CGAL::parameters::non_manifold()`
  */
  parameters::internal::Manifold_options manifold();

  /*!
  \ingroup PkgMesh3Parameters

  The function `parameters::non_manifold()` is used to drive the
  meshing algorithm for surfaces.
  It does not ensure that the surface of the output mesh is a manifold surface.
  The manifold property of the output mesh might nevertheless result from an appropriate
  choice of meshing criteria.
  \sa `CGAL::make_mesh_3()`
  \sa `CGAL::refine_mesh_3()`
  \sa `CGAL::parameters::manifold_with_boundary()`
  \sa `CGAL::parameters::manifold()`
  */
  parameters::internal::Manifold_options non_manifold();

  /*!
  \ingroup PkgMesh3Parameters

  The function `parameters::manifold_with_boundary()` is used to drive the
  meshing algorithm for surfaces.
  It ensures that the surface of the output mesh is a manifold surface which
  may have boundaries.
  The manifold property of the output mesh can be achieved only if the input surface
  is a manifold.
  Note that the meshing algorithm provably terminates only if the input
  sharp edges have been protected, using the
  feature protection (see \ref Mesh_3Protectionof0and1dimensionalExposed).

  \sa `CGAL::make_mesh_3()`
  \sa `CGAL::refine_mesh_3()`
  \sa `CGAL::parameters::non_manifold()`
  \sa `CGAL::parameters::manifold()`
  */
  parameters::internal::Manifold_options manifold_with_boundary();

/*!
\ingroup PkgMesh3Parameters

The function `parameters::exude()` allows the user to trigger a call to `exude_mesh_3()` in the
`make_mesh_3()` and `refine_mesh_3()` mesh generation functions.
It also allows the user to pass parameters
to the optimization function `exude_mesh_3()` through these mesh generation functions.

\cgalHeading{Parameters}

The parameters are named parameters. They are the same (i.e.\ they have the same
name and the same default values) as the parameters of `exude_mesh_3()`
function. See its manual page for further details.

\cgalHeading{Example}

\code{.cpp}
// Mesh generation with an exudation step
C3t3 c3t3 = make_mesh_3<c3t3>(domain,
                              criteria,
                              parameters::exude());

refine_mesh_3(c3t3,
              domain,
              criteria,
              parameters::exude(parameters::time_limit=10));
\endcode

\sa `CGAL::parameters::no_exude()`
\sa `CGAL::exude_mesh_3()`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`

*/
parameters::internal::Exude_options exude(
  double parameters::time_limit = 0,
  double parameters::sliver_bound = 0);

/*!
\ingroup PkgMesh3Parameters

The function `parameters::features()` provides a value of internal type `Features`
to specify if 0 and 1-dimensional features have to be taken into account.
The provided value is a default value that triggers the representation
of corners and curves in the mesh when the domain is a model
of `MeshDomainWithFeatures_3`.

Provides a `Features_options` value such that
0 and 1-dimensional input features are taken into account
if domain is a model of the refined concept `MeshDomainWithFeatures_3`.

\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`
\sa `CGAL::parameters::no_features()`

*/
parameters::internal::Features_options features();

/*!
\ingroup PkgMesh3Parameters

The function `parameters::lloyd()` allows the user to trigger a call of
`lloyd_optimize_mesh_3()` in the mesh generation functions
`make_mesh_3()` and `refine_mesh_3()`. It also allows the user to pass
parameters to the optimization function
`lloyd_optimize_mesh_3()` through these mesh generation functions.

\cgalHeading{Parameters}

The parameters are named parameters. They are the same (i.e.\ they have the same
name and the same default values) as the parameters of the `lloyd_optimize_mesh_3()`
function. See its manual page for further details.

\cgalHeading{Example}

\code{.cpp}
// Mesh generation with lloyd optimization step
C3t3 c3t3 = make_mesh_3<c3t3>(domain,
                              criteria,
                              parameters::lloyd());

refine_mesh_3(c3t3,
              domain,
              criteria,
              parameters::lloyd(parameters::time_limit=10));

\endcode

\sa `CGAL::parameters::no_lloyd()`
\sa `CGAL::lloyd_optimize_mesh_3()`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`

*/
parameters::internal::Lloyd_options lloyd(
double parameters::time_limit = 0,
std::size_t parameters::max_iteration_number = 0,
double parameters::convergence = 0.02,
double parameters::freeze_bound = 0.01,
bool parameters::do_freeze=true);

/*!
\ingroup PkgMesh3Parameters

The function `parameters::no_exude()` allows the user to tell the mesh generation functions
`make_mesh_3()` and `refine_mesh_3()` that no exudation must be done.

\cgalHeading{Example}

\code{.cpp}
// Mesh generation without exudation
C3t3 c3t3 = make_mesh_3<c3t3>(domain,
                              criteria,
                              parameters::no_exude());
\endcode

\sa `CGAL::parameters::exude()`
\sa `CGAL::exude_mesh_3()`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`

*/
parameters::internal::Exude_options no_exude();

/*!
\ingroup PkgMesh3Parameters

The function `parameters::no_features()` allows the user to prevent the handling
of 0 and 1-dimensional features. This is useful when the
domain is a model of `MeshDomainWithFeatures_3`
and the user does not want corners and curves
to be accurately represented
in the mesh.

Returns a `Features_options` value that prevents the mesh generator
to take into account 0 and 1-dimensional input features.

\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`
\sa `CGAL::parameters::features()`

*/
parameters::internal::Features_options no_features();

/*!
\ingroup PkgMesh3Parameters

The function `parameters::no_lloyd()` allows the user to tell the mesh generation functions
`make_mesh_3()` and `refine_mesh_3()` that no lloyd optimization must be done.

\cgalHeading{Example}

\code{.cpp}
// Mesh generation without lloyd optimization
C3t3 c3t3 = make_mesh_3<c3t3>(domain,
                              criteria,
                              parameters::no_lloyd());
\endcode

\sa `CGAL::parameters::lloyd()`
\sa `CGAL::lloyd_optimize_mesh_3()`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`

*/
parameters::internal::Lloyd_options no_lloyd();

/*!
\ingroup PkgMesh3Parameters

The function `parameters::no_odt()` allows the user to tell the mesh generation functions
`make_mesh_3()` and `refine_mesh_3()` that no odt optimization must be done.

\cgalHeading{Example}

\code{.cpp}
// Mesh generation without odt optimization
C3t3 c3t3 = make_mesh_3<c3t3>(domain,
                              criteria,
                              parameters::no_odt());
\endcode

\sa `CGAL::parameters::odt()`
\sa `CGAL::odt_optimize_mesh_3()`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`

*/
parameters::internal::Odt_options no_odt();

/*!
\ingroup PkgMesh3Parameters

The function `parameters::no_perturb()` allows the user to tell mesh generation global functions
`make_mesh_3()` and `refine_mesh_3()` that no perturbation must be done.

\cgalHeading{Example}

\code{.cpp}
// Mesh generation without perturbation
C3t3 c3t3 = make_mesh_3<c3t3>(domain,
                              criteria,
                              parameters::no_perturb());
\endcode

\sa `CGAL::parameters::perturb()`
\sa `CGAL::perturb_mesh_3()`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`

*/
parameters::internal::Perturb_options no_perturb();

/*!
\ingroup PkgMesh3Parameters

The function `parameters::odt()` allows the user to trigger a call to
`CGAL::odt_optimize_mesh_3()` in
`CGAL::make_mesh_3()` and `CGAL::refine_mesh_3()` mesh optimization functions. It also
allows the user to pass parameters to the optimization function
`odt_optimize_mesh_3()` through these mesh generation functions.

\cgalHeading{Parameters}

The parameters are named parameters. They are the same (i.e.\ they have the same
name and the same default values) as the parameters of `odt_optimize_mesh_3()`
function. See its manual page for further details.

\cgalHeading{Example}

\code{.cpp}
// Mesh generation with odt optimization step
C3t3 c3t3 = make_mesh_3<c3t3>(domain,
                              criteria,
                              parameters::odt());

refine_mesh_3(c3t3,
              domain,
              criteria,
              parameters::odt(parameters::time_limit=10));
\endcode

\sa `CGAL::parameters::no_odt()`
\sa `CGAL::odt_optimize_mesh_3()`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`

*/
parameters::internal::Odt_options odt(
double parameters::time_limit = 0,
std::size_t parameters::max_iteration_number = 0,
double parameters::convergence = 0.02,
double parameters::freeze_bound = 0.01,
bool parameters::do_freeze=true);

/*!
\ingroup PkgMesh3Parameters

The function `parameters::perturb()` allows the user to trigger a call to
`perturb_mesh_3()` in
`make_mesh_3()` and `refine_mesh_3()` mesh generation functions. It also
allows the user to pass parameters
to the optimization function `perturb_mesh_3()` through these mesh generation functions.

\cgalHeading{Parameters}

The parameters are named parameters. They are the same (i.e.\ they have the same
name and the same default values) as the parameters of `perturb_mesh_3()`
function. See its manual page for further details.

\cgalHeading{Example}

\code{.cpp}
// Mesh generation with a perturbation step
C3t3 c3t3 = make_mesh_3<c3t3>(domain,
                              criteria,
                              parameters::perturb());

refine_mesh_3(c3t3,
              domain,
              criteria,
              parameters::perturb(parameters::time_limit=10));

\endcode

\sa `CGAL::parameters::no_perturb()`
\sa `CGAL::perturb_mesh_3()`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`

*/
parameters::internal::Perturb_options perturb(
  double parameters::time_limit = 0,
  double parameters::sliver_bound = 0);

} /* namespace parameters */

} /* namespace CGAL */
