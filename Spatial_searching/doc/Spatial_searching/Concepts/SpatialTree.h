/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalConcept

The concept `SpatialTree` defines the requirements for a tree supporting 
both neighbor searching and approximate range searching. 

\cgalHasModel `CGAL::Kd_tree<Traits,Splitter,UseExtendedNode>`

*/

class SpatialTree {
public:

/// \name Types 
/// @{

/*!
Search traits. 
*/ 
typedef unspecified_type SearchTraits; 


/*!
Dimension tag.
*/
typedef unspecified_type D;
/*!
Point type. 
*/ 
typedef unspecified_type Point_d; 

/*!
Bidirectional const iterator with value type `Point_d` that allows 
to enumerate all points in the tree. 
*/ 
typedef unspecified_type iterator; 

/*!
Node handle. 
*/ 
typedef unspecified_type Node_handle; 

/*!
const node handle. 
*/ 
typedef unspecified_type Node_const_handle; 

/*!
const iterator with value type `const Point_d*`. 
*/ 
typedef unspecified_type Point_d_iterator; 

/*!
Splitter. 
*/ 
typedef unspecified_type Splitter; 

/*!
Distance. 
*/ 
typedef unspecified_type Distance; 

/// @} 

/// \name Creation 
/// @{

/*!

Constructs a tree on the elements from the sequence 
`[first,beyond)`. 

*/ 
template <class InputIterator> Tree(InputIterator first, InputIterator beyond, SearchTraits t); 

/// @} 

/// \name Operations 
/// @{

/*!
Reports the points that are approximately contained by `q`. The value type 
of `OutputIterator` must be `Point_d`. 
*/ 
template <class OutputIterator, class FuzzyQueryItem> 
OutputIterator search(OutputIterator it, FuzzyQueryItem q); 

/*!
Returns a const iterator to the first point in the tree. 
*/ 
iterator begin() const; 

/*!
Returns the appropriate past-the-end const iterator. 
*/ 
iterator end() const; 

/*!

Returns a handle to the root node of the tree. 

*/ 
Node_handle root(); 

/*!

Returns a const handle to the root node of the tree. 

*/ 
Node_const_handle root() const; 

/*!
Returns a const 
reference to the bounding box of the root node of the tree. 
*/ 
const Kd_tree_rectangle<SearchTraits::FT,D>& bounding_box() const; 

/*!

Returns the number of items that are stored in the tree. 

*/ 
unsigned int size() const; 

/// @}

}; /* end SpatialTree */
