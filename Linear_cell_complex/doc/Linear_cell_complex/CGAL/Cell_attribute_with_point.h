
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `Cell_attribute_with_point` represents an attribute containing a point and containing an information when `Info_` is different from `void`. This class can typically be used to associate a point to each 0-cell of a combinatorial or a generalized map.

\cgalModels `CellAttributeWithPoint`

\tparam LCC a model of the `LinearCellComplex` concept.
\tparam Info_ the type of the information contained in the attribute, `void` for no information. Equal to `void` by default.
\tparam Tag is `::Tag_true` to enable the storage of a \link Cell_attribute::Dart_handle `Dart_handle`\endlink of the associated cell, `::Tag_false` otherwise. Equal to `::Tag_true` by default.
\tparam OnMerge a functor called when two attributes are merged. Equal to `Null_functor` by default.
\tparam OnSplit a functor called when one attribute is split in two. Equal to `Null_functor` by default.

\sa `CGAL::Linear_cell_complex_min_items<d>`
\sa `CGAL::Linear_cell_complex_for_combinatorial_map<d,d2,LCCTraits,Items,Alloc>`
\sa `CGAL::Linear_cell_complex_for_generalized_map<d,d2,LCCTraits,Items,Alloc>`

*/
template< typename LCC, typename Info_, typename Tag, typename OnMerge, typename OnSplit >
class Cell_attribute_with_point : public CGAL::Cell_attribute<LCC,Info_,Tag,OnMerge,OnSplit> {
public:

/// \name Types
/// @{

/*!

*/
typedef LCC::Point Point;

/// @}

}; /* end Cell_attribute_with_point */
} /* end namespace CGAL */
