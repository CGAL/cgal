
namespace CGAL {

/*!
\ingroup PkgRangeSegmentTreesDSearchStructure

\brief A \f$ d\f$-dimensional range tree stores points and can be used to determine all
points that lie inside a given \f$ d\f$-dimensional interval.

\cgalHeading{Implementation}

The construction of a \f$ d\f$-dimensional range tree takes \f$ {O}(n\log n^{d-1})\f$ 
time. The points in 
the query window are reported in time \f$ {O}(k+{\log}^d n )\f$, where \f$ k\f$ 
is the number of reported points. 
The tree uses \f$ {O}(n\log n^{d-1})\f$ storage. 

*/
template< typename Data, typename Window, typename Traits >
class Range_tree_d {
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
container `Traits`. 
*/ 
typedef unspecified_type Traits; 

/// @} 

/// \name Creation 
/// @{

/*!

A range tree is 
constructed, such that the subtree of each vertex is of the 
same type prototype 
`sublayer_tree` is. 

We assume that the dimension of the tree is \f$ d\f$. This means, that 
`sublayer_tree` is a prototype of a \f$ d-1\f$-dimensional 
tree. All data items of the \f$ d\f$-dimensional range tree 
have container type `Data`. The query window of the 
tree has container type 
`Window`. `Traits` 
provides access to the corresponding data slots of container 
`Data` and `Window` for the \f$ d\f$-th 
dimension. The traits class `Traits` 
must at least provide all functions and type definitions 
as described in, for example, the reference page for 
`tree_point_traits`. 
The template class 
described there is fully generic and should fulfill the most 
requirements one can have. 
In order 
to generate a one-dimensional range tree instantiate `Tree_anchor<Data, Window> sublayer_tree` with the same template parameters (`Data` and 
`Window`) `Range_tree_d` is defined. In 
order to construct a two-dimensional range tree, create 
`Range_tree_d` with 
a one-dimensional `Range_tree_d` with the 
corresponding `Traits` class of the first dimension. 

\pre `Traits::Data==Data` and `Traits::Window==Window.` 
*/ 
Range_tree_d<Data, Window, Traits> 
r(Tree_base<Data, Window> sublayer_tree); 

/// @} 

/// \name Operations 
/// @{

/*!
The tree is constructed according to the data items in the 
sequence between the element pointed by iterator `first` and 
iterator `last`. The data items of the iterator must 
have type `Data`. 
\pre This function can only be called once. If it is the first call the tree is build and `true` is returned. Otherwise, nothing is done but a \cgal warning is given and `false` returned. 
*/ 
template<class ForwardIterator> 
bool make_tree(ForwardIterator first, ForwardIterator last); 

/*!

All elements that 
lay inside the \f$ d\f$-dimensional interval defined through 
`win` are placed in the sequence container of 
`OutputIterator`; the output iterator that points 
to the last location the function wrote to is returned. 
*/ 
template<class OutputIterator> 
OutputIterator window_query(Window win, OutputIterator result); 

/*!
The tree structure is checked. For each 
vertex the subtree is checked on being valid and it is checked 
whether the value of the `Key_type` of a vertex 
corresponds to the highest `Key_type` 
value of the left subtree. 
*/ 
bool is_valid(); 

protected:

/*!
returns `true`, if the 
data of `object` lies between the start and endpoint of 
interval `win`. Returns `false` otherwise. 
*/ 
bool is_inside(Window win, 
Data object); 

/*!
returns false. 
*/ 
bool is_anchor(); 

/// @}

}; /* end Range_tree_d */
} /* end namespace CGAL */
