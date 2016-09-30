
namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsClasses

The class `Dart` represents a <I>d</I>D dart in a combinatorial map.

\f$ \beta_i\f$ pointers are coded in a array of <I>d+1</I> \link Dart::Dart_handle `Dart_handle`\endlink
(because we describe also the \f$ \beta_0\f$ link). Attributes are
associated to each dart by \link Attribute_handle `Attribute_handle<i>`\endlink, one for each
non void <I>i</I>-attribute.

\cgalModels `::BasicDart`

\tparam d an integer for the dimension of the dart.

\tparam CMap must be a model of the `CombinatorialMap` concept.

\cgalHeading{Complexity}

\sa `CombinatorialMap`

*/
template< typename d, typename CMap >
class Dart {
public:

/// \name Types
/// @{

/*!

*/
typedef CMap::Dart_handle Dart_handle;

/*!

*/
typedef CMap::Dart_const_handle Dart_const_handle;

/*!

*/
template <unsigned int i>
using Attribute_handle = CMap::Attribute_handle<i>;

/*!

*/
template <unsigned int i>
using Attribute_const_handle = CMap::Attribute_const_handle<i>;

/// @}

/// \name Access Member Functions
/// @{

/*!
To simplify a future implementation, it is recommended to not use this function and to use \link CombinatorialMap::beta `cm.beta(dh,i)`\endlink instead.
Returns \f$ \beta_i\f$(`*this`).
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dimension</I>.
*/
Dart_handle beta(unsigned int i);

/*!
To simplify a future implementation, it is recommended to not use this function and to use \link CombinatorialMap::beta `cm.beta(dh,i)`\endlink instead.
Returns \f$ \beta_i\f$(`*this`) when the dart is const.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dimension</I>.
*/
Dart_const_handle beta(unsigned int i) const;

/*!
To simplify a future implementation, it is recommended to not use this function and to use \link CombinatorialMap::beta `cm.beta(dh,CGAL_BETAINV(i))`\endlink instead.
Returns \f$ \beta_i^{-1}\f$(`*this`).
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dimension</I>.
*/
Dart_handle beta_inv(unsigned int i);

/*!
To simplify a future implementation, it is recommended to not use this function and to use \link CombinatorialMap::beta `cm.beta(dh,CGAL_BETAINV(i))`\endlink instead.
Returns \f$ \beta_i^{-1}\f$(`*this`) when the dart is const.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dimension</I>.
*/
Dart_const_handle beta_inv(unsigned int i) const;

}; /* end Dart */
} /* end namespace CGAL */
