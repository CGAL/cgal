
namespace CGAL {

/*!
\ingroup PkgGeneralizedMapsClasses

The class `GDart` represents a <I>d</I>D dart in a generalized map.

\f$ \alpha_i\f$ pointers are coded in a array of <I>d+1</I> \ref GDart::Dart_handle "Dart_handle". Attributes are associated to each dart by \ref Attribute_handle "Attribute_handle<i>", one for each non void <I>i</I>-attribute.

\cgalModels `::GDart`

\tparam d an integer for the dimension of the dart.

\tparam GMap must be a model of the `GeneralizedMap` concept.

\cgalHeading{Complexity}

\sa `GeneralizedMap`

*/
template< typename d, typename GMap >
class GDart {
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

}; /* end GDart */
} /* end namespace CGAL */
