
namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsClasses

The class `Cell_attribute_with_id` represents an attribute containing (or not) an information, and having an id.

\cgalModels `CellAttribute`

\tparam Map a model of the `GenericMap` concept.

\tparam Info_ the type of the information contained in the attribute. Equal to `void` by default.

\tparam Tag is `::Tag_true` to enable the storage of a `Dart_handle` of the associated cell, `::Tag_false` otherwise. Equal to `::Tag_true` by default.

\tparam OnMerge the type of the functor called before two attributes are merged. Equal to `Null_functor` by default.

\tparam OnSplit the type of the functor called after one attribute is split in two. Equal to `Null_functor` by default.

\sa `CGAL::Combinatorial_map<d,Items,Alloc>`
\sa `CGAL::Generalized_map<d,Items,Alloc>`

*/
template< typename Map, typename Info_, typename Tag, typename OnMerge, typename OnSplit >
class Cell_attribute_with_id : public CGAL::Cell_attribute<Map,Info_,Tag,OnMerge,OnSplit> {
public:

/// Returns the index.
int id() const;

/// Returns the index by reference.
int& id();

}; /* end Cell_attribute_with_id */
} /* end namespace CGAL */
