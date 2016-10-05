
namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsClasses

The class `Cell_attribute` represents an attribute containing (or not) an information.

\cgalModels `CellAttribute`

\tparam Map a model of the `GenericcMap` concept.

\tparam Info_ the type of the information contained in the attribute.

\tparam Tag is `::Tag_true` to enable the storage of a `Dart_handle` of the associated cell, `::Tag_false` otherwise.

\tparam OnMerge the type of the functor called before two attributes are merged.

\tparam OnSplit the type of the functor called after one attribute is split in two.

By default, `OnMerge` and `OnSplit` are equal to `Null_functor`; `Tag` is equal to `::Tag_true`; and `Info_` is equal to `void`.

\sa `CGAL::Combinatorial_map<d,Items,Alloc>`
\sa `CGAL::Generalized_map<d,Items,Alloc>`

*/
template< typename Map, typename Info_, typename Tag, typename OnMerge, typename OnSplit >
class Cell_attribute {
public:

/// \name Types
/// @{

/*!

*/
typedef Info_ Info;

/*!

*/
typedef Tag Supports_cell_dart;

/*!

*/
typedef OnMerge On_merge;

/*!

*/
typedef OnSplit On_split;

/*!

*/
typedef CMap::Dart_handle Dart_handle;

/*!

*/
typedef CMap::Dart_const_handle Dart_const_handle;

/// @}

}; /* end Cell_attribute */
} /* end namespace CGAL */
