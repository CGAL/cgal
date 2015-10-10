namespace CGAL {

/*!
\ingroup AdvancedClasses

The class `Kd_tree_node` implements a node class for a `k-d` tree. 

A node is either a `Kd_tree_leaf_node` or a `Kd_tree_internal_node`.
A leaf node contains one or more points. An internal node contains 
a pointer to its lower child, a pointer to its
upper child, and a pointer to its separator. If extended nodes are 
used, it also contains the upper limit of the lower child's 
tight rectangle and the lower limit of the upper child's tight rectangle 
along the node's cutting dimension.


\tparam Traits must be a model of the concept `SearchTraits`, 
for example `Search_traits_2<Simple_cartesian<double> >`, or `Cartesian_d<double>`. 

\tparam Splitter must be a model of the concept `Splitter`.

*/
template< typename Traits, typename Splitter, typename UseExtendedNode >
class Kd_tree_node {
public:

/// \name Types 
/// @{
/*!
Point class.
*/
typedef Traits::Point_d Point_d;
/// @} 

/// \name Operations 
/// @{

/*!
Reports the points from the subtree of the node, that are approximately contained by `q`. 
*/ 
template <class OutputIterator, class FuzzyQueryItem> 
OutputIterator search(OutputIterator it, FuzzyQueryItem q) const; 

/*!
Reports any point from the subtree of the node, that is approximately contained by `q`. 
*/ 
template <class FuzzyQueryItem> 
boost::optional<Point_d> search_any_point(OutputIterator it, FuzzyQueryItem q) const; 

/*!
Reports all the points contained by the subtree of the node. 
*/ 
template <class OutputIterator> 
OutputIterator tree_items(OutputIterator it) const; 

/*!
Indicates whether a node is a leaf node. 
*/ 
bool is_leaf() const; 

/// @}

}; /* end Kd_tree_node */
/*!
\ingroup AdvancedClasses

*/
template < class TreeTraits, class Splitter, class UseExtendedNode > 
  class Kd_tree_leaf_node : public Kd_tree_node< TreeTraits, Splitter, UseExtendedNode >{
  public:

/// \name Types 
/// @{

/*!
const iterator over points. 
*/ 
typedef Kd_tree<Traits,Splitter,UseExtendedNode>::iterator iterator; 

/// @}

/// \name Operations 
/// @{

/*!
Returns the number of items stored in a leaf node. 
*/ 
unsigned int size() const; 

/*!
Returns a const iterator to the first item in a leaf node. 
*/ 
iterator begin() const; 

/*!
Returns the appropriate past-the-end const iterator. 
*/ 
iterator end() const; 

/// @}
}; /* end Kd_tree_leaf_node */

/*!
\ingroup AdvancedClasses

*/
template < class TreeTraits, class Splitter, class UseExtendedNode > 
  class Kd_tree_internal_node : public Kd_tree_node< TreeTraits, Splitter, UseExtendedNode >{
public:

/// \name Types 
/// @{

/*!
Number type. 
*/ 
typedef Traits::FT FT; 

/*!
Node handle. 
*/ 
typedef Kd_tree<Traits,Splitter,UseExtendedNode>::Node_handle Node_handle; 

/*!
const node handle. 
*/ 
typedef Kd_tree<Traits,Splitter,UseExtendedNode>::Node_const_handle Node_const_handle; 

/*!
Separator type. 
*/ 
typedef Splitter::Separator Separator; 

/// @}

/// \name Operations 
/// @{

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
Returns the cutting dimension. 
*/ 
int cutting_dimension() const;

/*!
Returns the cutting value. 
*/ 
int cutting_value() const;

/*!
Returns the upper limit of an extended node's 
lower child's tight rectangle 
along the node's cutting dimension. 
*/ 
FT low_value() const; 

/*!
Returns the lower limit of an extended node's
upper child's tight rectangle 
along the node's cutting dimension. 
*/ 
FT high_value() const; 

/// @}
}; /* end Kd_tree_leaf_node */

} /* end namespace CGAL */
