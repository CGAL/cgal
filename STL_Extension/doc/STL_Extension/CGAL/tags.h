
namespace CGAL {

/*!
\ingroup PkgStlExtension

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
\ingroup PkgStlExtension

The typedef `Tag_false` is `Boolean_tag<false>`. 
It is used to indicate, for example, 
that a certain feature is not available in a class. 

\sa `CGAL::Boolean_tag<bool value>` 
\sa `CGAL::Tag_true` 
*/
class Tag_false {
public:

/// \name Definition 
/// @{ 
/*! 
is `false` 
*/ 
static const bool value; 

/// @} 

}; /* end Tag_false */


/*!
\ingroup PkgStlExtension

The typedef `Tag_true` is `Boolean_tag<true>`. 
It is used to indicate, for example, 
that a certain feature is available in a class. 

\sa `CGAL::Boolean_tag<bool value>` 
\sa `CGAL::Tag_false` 


*/
class Tag_true {
public:

/// \name Definition 
/// @{ 
/*! 
is `true` 
*/ 
static const bool value; 

/// @} 

}; /* end Tag_true */

/*!
\ingroup PkgStlExtension

Class indicating the absence of a functor.
\models ::DefaultConstructible 

\sa `AlgebraicStructureTraits` 
\sa `RealEmbeddableTraits` 
*/
struct Null_functor {
  
};

/*!
\ingroup PkgStlExtension

General tag indicating that non of any other possible tags is valid. 

\models ::DefaultConstructible 

\sa `AlgebraicStructureTraits` 

*/
class Null_tag {
public:
}; /* end Null_tag */

} /* end namespace CGAL */
