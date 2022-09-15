namespace CGAL {

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
  unspecified_type manifold();

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
  unspecified_type non_manifold();

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
  unspecified_type manifold_with_boundary();

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
unspecified_type exude(
  double parameters::time_limit = 0,
  double parameters::sliver_bound = 0);

/*!
\ingroup PkgMesh3Parameters

The function `parameters::features()` can be used to specify
that 0 and 1-dimensional features have to be taken into account.
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
unspecified_type features();

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
unspecified_type lloyd(
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
unspecified_type no_exude();

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
unspecified_type no_features();

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
unspecified_type no_lloyd();

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
unspecified_type no_odt();

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
unspecified_type no_perturb();

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
unspecified_type odt(
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
unspecified_type perturb(
  double parameters::time_limit = 0,
  double parameters::sliver_bound = 0);

} /* namespace parameters */

} /* namespace CGAL */
