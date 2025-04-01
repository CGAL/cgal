
namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities

Depending on `bool value` the class `Boolean_tag` indicates that
something is `true` or `false` respectively.

\sa `CGAL::Tag_true`
\sa `CGAL::Tag_false`

*/
template< typename bool value >
struct Boolean_tag {

/// \name Constants
/// @{
/*!

*/
static const bool value;

/// @}

}; /* end Boolean_tag */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities

The typedef `Tag_false` is `Boolean_tag<false>`.
It is used to indicate, for example,
that a certain feature is not available in a class.

\sa `CGAL::Boolean_tag<bool value>`
\sa `CGAL::Tag_true`
*/
typedef CGAL::Boolean_tag<false> Tag_false;


/*!
\ingroup PkgSTLExtensionUtilities

The typedef `Tag_true` is `Boolean_tag<true>`.
It is used to indicate, for example,
that a certain feature is available in a class.

\sa `CGAL::Boolean_tag<bool value>`
\sa `CGAL::Tag_false`
*/
typedef CGAL::Boolean_tag<true> Tag_true;


/*!
\ingroup PkgSTLExtensionUtilities

Class indicating the absence of a functor.
\cgalModels{DefaultConstructible}

\sa `AlgebraicStructureTraits`
\sa `RealEmbeddableTraits`
*/
struct Null_functor {

};

/*!
\ingroup PkgSTLExtensionUtilities
Tag used to disable concurrency.
For example, it may be used by a user to request the sequential version of an algorithm.
*/
struct Sequential_tag {};

/*!
\ingroup PkgSTLExtensionUtilities
Tag used to enable concurrency.
For example, it may be used by a user to request the parallel version of an algorithm.
*/
struct Parallel_tag {};

/*!
\ingroup PkgSTLExtensionUtilities
This tag is a convenience typedef to `Parallel_tag` if the third party library \ref thirdpartyTBB
has been found and linked, and to `Sequential_tag` otherwise.
*/
struct Parallel_if_available_tag {};

/*!
\ingroup PkgSTLExtensionUtilities

General tag indicating that non of any other possible tags is valid.

\cgalModels{DefaultConstructible}

\sa `AlgebraicStructureTraits`

*/
struct Null_tag {
}; /* end Null_tag */

} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities

The class `Manifold_tag` is a tag class used to monitor the
surface meshing algorithm. When instantiated with the tag
`Manifold_tag` the function template
`make_surface_mesh()`
ensures that the output mesh is a manifold surface
without boundary.

\sa `make_surface_mesh()`
\sa `Manifold_with_boundary_tag`
\sa `Non_manifold_tag`

*/

struct Manifold_tag {

}; /* end Manifold_tag */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities

The class `Manifold_with_boundary_tag` is a tag class used to monitor the
surface meshing algorithm. When instantiated with the tag
`Manifold_with_boundary_tag`, the function template
`make_surface_mesh()`
ensures that the output mesh is a manifold surface
but it may have boundaries.

\sa `make_surface_mesh()`
\sa `Manifold_tag`
\sa `Non_manifold_tag`

*/

struct Manifold_with_boundary_tag {

}; /* end Manifold_with_boundary_tag */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities

The class `Non_manifold_tag` is a tag class used to monitor the
surface meshing algorithm. When instantiated with the tag
`Non_manifold_tag` the function template
`make_surface_mesh()`
does not ensure that the output mesh is a manifold surface.
The manifold property of output mesh
may nevertheless result from the choice of
appropriate meshing criteria.

\sa `make_surface_mesh()`
\sa `Manifold_tag`
\sa `Manifold_with_boundary_tag`

*/

struct Non_manifold_tag {

}; /* end Non_manifold_tag */
} /* end namespace CGAL */
