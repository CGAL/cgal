
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `GMap_linear_cell_complex` represents a linear cell complex in dimension `d`, in an ambient space of dimension `d2`, using a generalized map as underlying combinatorial data-structure. This is a model of the concept of `LinearCellComplex` and a model of the concept `GeneralizedMap`.

\cgalModels `LinearCellComplex`
\cgalModels `GeneralizedMap`

\tparam d an integer for the dimension of the generalized map,
\tparam d2 an integer for the dimension of the ambient space,
\tparam LCCTraits must be a model of the `LinearCellComplexTraits` concept, satisfying \link LinearCellComplexTraits::ambient_dimension `LCCTraits::ambient_dimension`\endlink`==d2`,
\tparam Items must be a model of the `LinearCellComplexItems` concept,
\tparam Alloc has to match the standard allocator requirements.

There are four default template arguments: `d2` is equal to `d`, `LCCTraits` is equal to `CGAL::Linear_cell_complex_traits<d2>`, `Items` is equal to `CGAL::Linear_cell_complex_min_items<d>` and `Alloc` is `CGAL_ALLOCATOR(int)`.

\cgalAdvancedBegin
Note that there is an additional, and undocumented, template parameter `GMap` for `Linear_cell_complex<d,d2,LCCTraits,Items,Alloc,GMap>` allowing to inherit from any model of the `GeneralizedMap` concept. \cgalAdvancedEnd

\sa `CGAL::Generalized_map<d,Items,Alloc>`
\sa `CGAL::Linear_cell_complex_traits<d,K>`
\sa `CGAL::GMap_linear_cell_complex_min_items<d>`

*/

template< typename d, typename d2, typename LCCTraits, typename Items, typename Alloc >
class GMap_linear_cell_complex  : public Generalized_map<d,Items,Alloc>
{
public:

/// \name Constants
/// @{

/*!
Ambient dimension, must be > 1.
*/
static unsigned int ambient_dimension = d2;

/// @}

/// \name Types
/// @{

/*!

*/
typedef GMap_linear_cell_complex<d,d2,LCCTraits,Items,Alloc> Self;

/*!
The type of dart, must satisfy \link Dart::dimension `Dart::dimension`\endlink`==d`.
*/
typedef Items::Dart_wrapper<Self>::Dart Dart;
  
}; /* end GMap_linear_cell_complex */
  
} /* end namespace CGAL */
