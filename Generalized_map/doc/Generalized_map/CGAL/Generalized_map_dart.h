
namespace CGAL {

/*!
\ingroup PkgGeneralizedMapsClasses

The class `Generalized_map_dart` represents a <I>d</I>D dart in a generalized map.

\f$ \alpha_i\f$ pointers are coded in a array of <I>d+1</I> \link Generalized_map_dart::Dart_handle `Dart_handle`\endlink. Attributes are associated to each dart by \link Attribute_handle `Attribute_handle<i>`\endlink, one for each non void <I>i</I>-attribute.

\cgalModels `Dart`

\tparam d the dimension of the dart.

\tparam GMap a model of the `GeneralizedMap` concept.

\sa `GeneralizedMap`

*/
template< typename d, typename GMap >
class Generalized_map_dart {
public:

/// \name Types
/// @{

/*!

*/
typedef GMap::Dart_handle Dart_handle;

/*!

*/
typedef GMap::Dart_const_handle Dart_const_handle;

/*!

*/
template <unsigned int i>
using Attribute_handle = GMap::Attribute_handle<i>;

/*!

*/
template <unsigned int i>
using Attribute_const_handle = GMap::Attribute_const_handle<i>;

/// @}

}; /* end Generalized_map_dart */
} /* end namespace CGAL */
