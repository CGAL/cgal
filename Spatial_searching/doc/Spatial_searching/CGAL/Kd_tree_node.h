namespace CGAL {

/*!
\ingroup AdvancedClasses

The class `Kd_tree_node` implements a node class for a `k-d` tree. 

\cgalAdvanced A node is either a leaf node, an internal node or an
extended internal node.  A leaf node contains one or more points. An
internal node contains a pointer to its lower child, a pointer to its
upper child, and a pointer to its separator.  An extended internal
node is an internal node containing the lower and upper limit of an
extended node's rectangle along the node's cutting dimension.

### Parameters ###

Expects for the template argument a model of the concept `SearchTraits`, 
for example `CGAL::Search_traits_2<CGAL::Cartesian<double> >`, or `CGAL::Cartesian_d<double>`. 

*/
template< typename Traits, typename Splitter, typename UseExtendedNode >
class Kd_tree_node {
public:

/// \name Types 
/// @{

/*! 
Denotes type of node. 
*/ 

/*! 
Number type. 
*/ 
typedef Traits::FT FT; 

/*! 
Point type. 
*/ 
typedef Traits::Point_d Point_d; 

/*! 
Separator type. 
*/ 
typedef Splitter::Separator Separator; 

/*! 
const iterator over points. 
*/ 
typedef Kd_tree<Traits,Splitter,UseExtendedNode>::Point_d_iterator Point_d_iterator; 

/*! 
Node handle. 
*/ 
typedef Kd_tree<Traits,Splitter,UseExtendedNode>::Node_handle Node_handle; 

/*! 
const node handle. 
*/ 
typedef Kd_tree<Traits,Splitter,UseExtendedNode>::Node_const_handle Node_const_handle; 

/// @} 

/// \name Operations 
/// @{

/*! 
Reports the points from the subtree of the node, that are approximately contained by q. 
*/ 
template <class OutputIterator, class FuzzyQueryItem> 
OutputIterator search(OutputIterator it, FuzzyQueryItem q) const; 

/*! 
Reports all the points contained by the subtree of the node. 
*/ 
template <class OutputIterator> 
OutputIterator tree_items(OutputIterator it) const; 

/*! 
Indicates whether a node is a leaf node. 
*/ 
bool is_leaf() const; 

/*! 
Returns the number of items stored in a leaf node. 
*/ 
unsigned int size() const; 

/*! 
Returns a const iterator to the first item in a leaf node. 
*/ 
Point_d_iterator begin() const; 

/*! 
Returns the appropriate past-the-end const iterator. 
*/ 
Point_d_iterator end() const; 

/*! 
Returns a handle to the lower child of an internal node. 
*/ 
Node_handle lower(); 

/*! 
Returns a handle to the upper child of an internal node. 
*/ 
Node_handle upper(); 

/*! 
Returns a const handle to the lower child of an internal node. 
*/ 
Node_const_handle lower() const; 

/*! 
Returns a const handle to the upper child of an internal node. 
*/ 
Node_const_handle upper() const; 

/*! 
Returns a reference to the separator. 
*/ 
Separator& separator(); 

/*! 
Returns the lower limit of an extended node's rectangle 
along the node's cutting dimension. 
*/ 
FT low_value() const; 

/*! 
Returns the upper limit of an extended node's rectangle 
along the node's cutting dimension. 
*/ 
FT high_value() const; 

/// @}

}; /* end Kd_tree_node */
} /* end namespace CGAL */
