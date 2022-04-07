
namespace CGAL {

/*!
\ingroup PkgSearchStructuresSearchStructure

\brief  A \f$ d\f$-dimensional segment tree stores  \f$ d\f$-dimensional intervals and can be used to find all intervals that enclose, partially overlap, or contain a query interval, which may be a point.

\cgalHeading{Implementation}

A \f$ d\f$-dimensional segment tree is constructed in \f$ {O}(n\log n^d)\f$ time.
An inverse range query is performed in time \f$ {O}(k+{\log}^d n )\f$, where \f$ k\f$
is the number of reported intervals.
The tree uses \f$ {O}(n\log n^d)\f$ storage.

*/
template< typename Data, typename Window, typename Traits >
class Segment_tree_d {
public:

/// \name Types
/// @{

/*!
container `Data`.
*/
typedef unspecified_type Data;

/*!
container `Window`.
*/
typedef unspecified_type Window;

/*!
class `Traits`.
*/
typedef unspecified_type Traits;

/// @}

/// \name Creation
/// @{

/*!
A segment tree is defined, such that the subtree of each
vertex is of the same type prototype
`sublayer_tree` is.

We assume that the dimension of the tree is \f$ d\f$. This means, that
`sublayer_tree` is a prototype of a \f$ d-1\f$-dimensional
tree. All data items of the \f$ d\f$-dimensional segment tree
have container type `Data`. The query window of the
tree has container type
`Window`. `Traits`
provides access to the corresponding data slots of container
`Data` and `Window` for the \f$ d\f$-th
dimension. The traits class `Traits`
must at least provide all functions and type definitions
described, for example, in the reference page for
`tree_point_traits`.
The template class
described there is fully generic and should fulfill the most
requirements one can have.
In order
to generate a one-dimensional segment tree instantiate `Tree_anchor<Data, Window> sublayer_tree` with the same template parameters `Data` and
`Window` `Segment_tree_d` is defined. In
order to construct a two-dimensional segment tree, create
`Segment_tree_d` with
a one-dimensional `Segment_tree_d` with the
corresponding `Traits` of the first dimension.

\pre `Traits::Data==Data` and `Traits::Window==Window.`
*/
Segment_tree_d<Data, Window,
Traits> s(Tree_base<Data, Window> sublayer_tree);

/// @}

/// \name Operations
/// @{

/*!
The tree is constructed according to the data items in the
sequence between the element pointed by iterator `first` and
iterator `last`.
\pre This function can only be called once. If it is the first call the tree is build and `true` is returned. Otherwise, nothing is done but a `CGAL warning` is given and `false` returned.
*/
bool make_tree(In_it first, In_it last);

/*!
`win`\f$ =[a_1,b_1),\ldots, [a_d,b_d)\f$, \f$ a_i,b_i\in T_i\f$, \f$ 1\le i\le d\f$. All elements that
intersect the associated \f$ d\f$-dimensional interval of
`win` are placed in the
associated sequence container of `OutputIterator` and
returns an
output iterator that points
to the last location the function wrote to.
In order to perform an inverse range query, a range query of
\f$ \epsilon\f$ width has to be performed.

*/
OutputIterator window_query(Window win, OutputIterator result);

/*!
All elements that
enclose the associated \f$ d\f$-dimensional interval of
`win` are placed in the
associated sequence container of `OutputIterator` and returns an output iterator that points
to the last location the function wrote to.
*/
OutputIterator enclosing_query(Window win, OutputIterator result);

/*!
The tree structure is checked. For each
vertex either the
sublayer tree is a tree anchor, or it stores a (possibly empty)
list of data items. In the first case, the sublayer tree of the
vertex is checked on being valid. In the second case, each data
item is checked weather it contains the associated interval of
the vertex and does not contain the associated interval of the
parent vertex or not. `true` is returned if the tree structure is valid,
`false` otherwise.
*/
bool is_valid();

protected:

/*!
returns `true`, if the
interval of `object` is contained in the
interval of `win`, `false` otherwise.
*/
bool is_inside(Window win, Data object);

/*!
returns false.
*/
bool is_anchor();

/// @}

}; /* end Segment_tree_d */
} /* end namespace CGAL */
