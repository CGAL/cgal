
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `Linear_cell_complex_for_generalized_map` represents a linear cell complex in dimension `d`, in an ambient space of dimension `d2`, using a generalized map as underlying combinatorial data-structure.

\cgalModels `LinearCellComplex`
\cgalModels `GeneralizedMap`

\tparam d the dimension of the generalized map.
\tparam d2 the dimension of the ambient space. Equal to  `d` by default.
\tparam LCCTraits a model of the `LinearCellComplexTraits` concept, satisfying \link LinearCellComplexTraits::ambient_dimension `LCCTraits::ambient_dimension`\endlink`==d2`. Equal to`CGAL::Linear_cell_complex_traits<d2>` by default.
\tparam Items a model of the `LinearCellComplexItems` concept. Equal to `CGAL::Linear_cell_complex_min_items<d>` by default.
\tparam Alloc has to match the standard allocator requirements. Equal to `CGAL_ALLOCATOR(int)` by default.

\cgalAdvancedBegin
Note that there is an additional, and undocumented, template parameter `GMap` for `Linear_cell_complex_for_generalized_map<d,d2,LCCTraits,Items,Alloc,GMap>` allowing to inherit from any model of the `GeneralizedMap` concept. \cgalAdvancedEnd

\sa `CGAL::Generalized_map<d,Items,Alloc>`
\sa `CGAL::Linear_cell_complex_traits<d,K>`
\sa `CGAL::Linear_cell_complex_min_items<d>`

*/

template< typename d, typename d2, typename LCCTraits, typename Items, typename Alloc >
class Linear_cell_complex_for_generalized_map  : public Generalized_map<d,Items,Alloc>
{
public:

/// \name Constants
/// @{

/*!
Ambient dimension, must be > 1.
*/
static unsigned int ambient_dimension = d2;

/// @}

}; /* end Linear_cell_complex_for_generalized_map */

} /* end namespace CGAL */
