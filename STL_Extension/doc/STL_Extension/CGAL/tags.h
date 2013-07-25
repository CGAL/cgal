
namespace CGAL {

/*!
\ingroup PkgStlExtensionUtilities

Depending on `bool value` the class `Boolean_tag` indicates that 
something is `true` or `false` respectively. 

\sa `CGAL::Tag_true` 
\sa `CGAL::Tag_false` 

*/
template< typename bool value >
class Boolean_tag {
public:

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
\ingroup PkgStlExtensionUtilities

The typedef `Tag_false` is `Boolean_tag<false>`. 
It is used to indicate, for example, 
that a certain feature is not available in a class. 

\sa `CGAL::Boolean_tag<bool value>` 
\sa `CGAL::Tag_true` 
*/
typedef CGAL::Boolean_tag<false> Tag_false;


/*!
\ingroup PkgStlExtensionUtilities

The typedef `Tag_true` is `Boolean_tag<true>`. 
It is used to indicate, for example, 
that a certain feature is available in a class. 

\sa `CGAL::Boolean_tag<bool value>` 
\sa `CGAL::Tag_false` 
*/
typedef CGAL::Boolean_tag<true> Tag_true;


/*!
\ingroup PkgStlExtensionUtilities

Class indicating the absence of a functor.
\cgalModels `DefaultConstructible`

\sa `AlgebraicStructureTraits` 
\sa `RealEmbeddableTraits` 
*/
struct Null_functor {
  
};

/*!
\ingroup PkgStlExtensionUtilities
Tag used to enable/disable concurrency.
For example, it may be used by a user to request the sequential version of an algorithm.
*/
struct Sequential_tag {};

/*!
\ingroup PkgStlExtensionUtilities
Tag used to enable/disable concurrency.
For example, it may be used by a user to request the parallel version of an algorithm.
*/
struct Parallel_tag {};

/*!
\ingroup PkgStlExtensionUtilities

General tag indicating that non of any other possible tags is valid. 

\cgalModels `DefaultConstructible`

\sa `AlgebraicStructureTraits` 

*/
class Null_tag {
public:
}; /* end Null_tag */

} /* end namespace CGAL */
