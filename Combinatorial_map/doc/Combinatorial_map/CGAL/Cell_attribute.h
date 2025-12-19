
namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsClasses

The class `Cell_attribute` represents an attribute containing (or not) an information.

\cgalModels{CellAttribute}

\tparam Map a model of the `GenericMap` concept.

\tparam Info_ the type of the information contained in the attribute. Equal to `void` by default.

\tparam Tag is `::Tag_true` to enable the storage of a `Dart_descriptor` of the associated cell, `::Tag_false` otherwise. Equal to `::Tag_true` by default.

\tparam OnMerge the type of the functor called before two attributes are merged. Equal to `Null_functor` by default.

\tparam OnSplit the type of the functor called after one attribute is split in two. Equal to `Null_functor` by default.

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
typedef CMap::Dart_descriptor Dart_descriptor;

/*!

*/
typedef CMap::Dart_const_descriptor Dart_const_descriptor;

/// @}

}; /* end Cell_attribute */
} /* end namespace CGAL */
