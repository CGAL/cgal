
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `Linear_cell_complex_for_combinatorial_map` represents a linear cell complex in dimension `d`, in an ambient space of dimension `d2`, using a combinatorial map as underlying combinatorial data-structure.

\cgalModels `LinearCellComplex`
\cgalModels `CombinatorialMap`

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

\deprecated Before CGAL 4.9, this class was named `%Linear_cell_complex`. This old name still exist for backward compatibility.

\deprecated Before CGAL 4.9, `Items` had to define the type of dart used. This is now deprecated, the `Dart` type is no more defined in the item class, but replaced by the `Dart_info` type. See deprecated note in the `Linear_cell_complex_min_items` class. `CGAL_CMAP_DART_DEPRECATED` can be defined to keep the old behavior.

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
