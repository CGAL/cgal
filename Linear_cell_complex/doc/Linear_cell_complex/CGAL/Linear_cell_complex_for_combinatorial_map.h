
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `Linear_cell_complex_for_combinatorial_map` represents a linear cell complex in dimension `d`, in an ambient space of dimension `d2`, using a combinatorial map as underlying combinatorial data-structure.
Like for `Combinatorial_map`, two versions exist: one where Darts and non void attributes are stored in memory using `Compact_container`, using `Alloc` as allocator, and use handles as descriptors; a second one where Darts and non void attributes are stored in an internal std::vector like data-structure, and use indices as descriptors. The choice between the two versions is done through the item class.

\cgalModels{LinearCellComplex,CombinatorialMap}

\tparam d the dimension of the combinatorial map.
\tparam d2 the dimension of the ambient space. Equal to `d` by default.
\tparam LCCTraits be a model of the `LinearCellComplexTraits` concept, satisfying \link LinearCellComplexTraits::ambient_dimension `LCCTraits::ambient_dimension`\endlink`==d2`. Equal to `CGAL::Linear_cell_complex_traits<d2>` by default.
\tparam Items a model of the `LinearCellComplexItems` concept. Equal to `CGAL::Linear_cell_complex_min_items` by default.
\tparam Alloc has to match the standard allocator requirements. Equal to `CGAL_ALLOCATOR(int)` by default.

\cgalAdvancedBegin
Note that there is an additional, and undocumented, template parameter `CMap` for `Linear_cell_complex_for_combinatorial_map<d,d2,LCCTraits,Items,Alloc,CMap>` allowing to inherit from any model of the `CombinatorialMap` concept. \cgalAdvancedEnd

\sa `CGAL::Combinatorial_map<d,Items,Alloc>`
\sa `CGAL::Linear_cell_complex_traits<d,K>`
\sa `CGAL::Linear_cell_complex_min_items<d>`

*/

template< typename d, typename d2, typename LCCTraits, typename Items, typename Alloc >
class Linear_cell_complex_for_combinatorial_map : public Combinatorial_map<d,Items,Alloc>
 {
public:

/// \name Constants
/// @{

/*!
Ambient dimension, must be > 1.
*/
static const unsigned int ambient_dimension = d2;

/// @}

}; /* end Linear_cell_complex_for_combinatorial_map */

} /* end namespace CGAL */
