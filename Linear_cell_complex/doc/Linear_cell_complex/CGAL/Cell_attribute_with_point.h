
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

The class `Cell_attribute_with_point` represents an attribute containing a point and
containing an information when `Info_` is different from `void`.
This class can typically be used to associate a point to each 0-cell
of a combinatorial map.

\cgalModels `CellAttributeWithPoint`

\tparam LCC must be an instantiation of `Linear_cell_complex` class,
\tparam Info_ is the type of the information contained in the attribute, `void` for no information,

\tparam Tag is `::Tag_true` to enable the storage of a \ref Cell_attribute_with_point::Dart_handle "Dart_handle" of the associated cell, `::Tag_false` otherwise,
\tparam OnMerge is a functor called when two attributes are merged,
\tparam OnSplit is a functor called when one attribute is split in two.

By default, `OnMerge` and `OnSplit` are equal to
`Null_functor`; `Tag` is equal to
`::Tag_true`; and `Info_` is equal to `void`.

\sa `CGAL::Linear_cell_complex<d,d2,LCCTraits,Items,Alloc>`
\sa `CGAL::Linear_cell_complex_min_items<d>`
\sa `CGAL::Cell_attribute<CMap,Info_,Tag,OnMerge,OnSplit>`

*/
template< typename LCC, typename Info_, typename Tag, typename OnMerge, typename OnSplit >
class Cell_attribute_with_point : public CGAL::Cell_attribute<CMap,Info_,Tag,OnMerge,OnSplit> {
public:

/// \name Types
/// @{

/*!

*/
typedef LCC::Point Point;

/*!

*/
typedef LCC::Dart_handle Dart_handle;

/*!

*/
typedef LCC::Dart_const_handle Dart_const_handle;

/// @}

}; /* end Cell_attribute_with_point */
} /* end namespace CGAL */
