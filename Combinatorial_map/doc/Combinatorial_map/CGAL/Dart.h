
namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsClasses

The class `Dart` represents a <I>d</I>D dart.

\f$ \beta_i\f$ pointers are coded in a array of <I>d+1</I> \ref Dart::Dart_handle "Dart_handle"
(because we describe also the \f$ \beta_0\f$ link). Attributes are
associated to each dart by \ref Attribute_handle "Attribute_handle<i>", one for each
non void <I>i</I>-attribute.

\cgalModels `::Dart`

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

}; /* end Dart */
} /* end namespace CGAL */
