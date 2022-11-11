
namespace CGAL {

/*!
\ingroup BGLGraphExternalIndices

The class `Linear_cell_complex_for_bgl_combinatorial_map_helper` defines a `CGAL::Linear_cell_complex_for_combinatorial_map` as inner type, named `type`, having `CGAL::Linear_cell_complex_bgl_min_items` as items class. With this item class, no information are associated with darts, darts have ids and 0- and 2-attributes are enabled and have ids.

\tparam d the dimension of the combinatorial map.
\tparam d2 the dimension of the ambient space. Equal to `d` by default.
\tparam LCCTraits be a model of the `LinearCellComplexTraits` concept, satisfying \link LinearCellComplexTraits::ambient_dimension `LCCTraits::ambient_dimension`\endlink`==d2`. Equal to `CGAL::Linear_cell_complex_traits<d2>` by default.
\tparam Alloc has to match the standard allocator requirements. Equal to `CGAL_ALLOCATOR(int)` by default.

\sa `CGAL::Linear_cell_complex_bgl_min_item`
\sa `CGAL::Linear_cell_complex_for_combinatorial_map<d,d2,Traits,Items,Alloc>`

*/
template< typename d, typename d2, typename LCCTraits, typename Alloc >
struct Linear_cell_complex_for_bgl_combinatorial_map_helper {
  /// Type of the Linear_cell_complex_for_combinatorial_map.
  typedef CGAL::Linear_cell_complex_for_combinatorial_map
            <d, d2, LCCTraits, CGAL::Linear_cell_complex_bgl_min_items,
             Alloc> type;

}; /* end Linear_cell_complex_for_bgl_combinatorial_map_helper */

} /* end namespace CGAL */
