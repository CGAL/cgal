
namespace CGAL {

/*!
\ingroup PkgRangeSegmentTreesD

`Tree_anchor` is also derived from `Tree_base`. Therefore, it provides
the same methods as `Range_tree_d` and `Segment_tree_d`, but does
nothing; it can be used as a recursion anchor for those
classes. Therefore, instantiate `Sublayer_type` of `Range_tree_d`
(`Segment_tree_d` respectively) with `Tree_anchor` and the container
classes for the data items (`Data` and `Window`).

\cgalHeading{Example}

The following figures show a number of rectangles and a 2-dimensional 
segment tree built on them. 

\image html segment_ex2.png "Two dimensional interval data and the corresponding segment tree."
\image latex segment_ex2.png "Two dimensional interval data and the corresponding segment tree."
*/
template< typename Data, typename Window >
class Tree_anchor : public Tree_base {
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

/// @} 

/// \name Creation 
/// @{

/*!

*/ 
Tree_anchor<Data, Window> a(); 

/// @} 

/// \name Operations 
/// @{

/*!

*/ 
template<class OutputIterator> 
OutputIterator window_query(Window win, OutputIterator result); 

/*!

*/ 
template<class OutputIterator> 
OutputIterator enclosing_query(Window win, OutputIterator result); 

/*!
returns true; 
*/ 
bool is_valid(); 

protected:

/*!
returns true. 
*/ 
bool is_inside(Window win, 
Data object); 

/*!
returns true. 
*/ 
bool is_anchor(); 

/// @}

}; /* end Tree_anchor */
} /* end namespace CGAL */
