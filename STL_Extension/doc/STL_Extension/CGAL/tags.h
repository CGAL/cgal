
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
