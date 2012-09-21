namespace CGAL {

/*!
\ingroup SearchClasses

The class `Kd_tree` defines a \f$ k\f$-\f$ d\f$ tree. 

### Parameters ###

Expects for the first template argument a model of the concept 
`SearchTraits`, for example `CGAL::Search_traits_2<CGAL::Cartesian<double> >`. 

Expects for the second template argument a model for the concept `Splitter`. 
It defaults to `Sliding_midpoint<Traits>`. 

Expects for the third template argument `CGAL::Tag_true`, if the 
tree shall be built with extended nodes, and `CGAL::Tag_false` otherwise. 

\sa Tree
\sa `CGAL::Kd_tree_node<Traits>` 
\sa `CGAL::Search_traits_2<Kernel>` 
\sa `CGAL::Search_traits_3<Kernel>` 
\sa `CGAL::Search_traits<FT_,Point,CartesianIterator,ConstructCartesianIterator>` 

*/
template< typename Traits, typename Splitter, typename UseExtendedNode >
class Kd_tree {
public:

/// \name Types 
/// @{

/*! 
Point class. 
*/ 
typedef Traits::Point_d Point_d; 

/*! 
Number type. 
*/ 
typedef Traits::FT FT; 

/*! 
Splitter type. 
*/ 
typedef Hidden_type Splitter; 

/*! 
Bidirectional const iterator with value type `Point_d` that allows 
to enumerate all points in the tree. 
*/ 
typedef Hidden_type iterator; 

/*! 
  \advanced A handle with value type `Kd_tree_node<Traits,Splitter>`. 
*/ 
typedef Hidden_type Node_handle; 

/*! 
  \advanced A const handle with value type `Kd_tree_node<Traits,Splitter>`. 
*/ 
typedef Hidden_type Node_const_handle; 

/*! 
  \advanced Random access const iterator with value type `const Point_d*`. 
*/ 
typedef Hidden_type Point_d_iterator; 

/*! 
  \advanced A type that counts the number of elements in a \f$ k\f$-\f$ d\f$ tree. 
*/ 
typedef Hidden_type size_type; 

/// @} 

/// \name Creation 
/// @{

/*! 
Constructs an empty \f$ k\f$-\f$ d\f$ tree. 
*/ 
Kd_tree(Splitter s=Splitter(),Traits t=Traits()); 

/*! 

Constructs a \f$ k\f$-\f$ d\f$ tree on the elements from the sequence 
`[first, beyond)` using the splitting rule implemented by `s`. 
The value type of the `InputIterator` must be `Point_d`. 

*/ 
template <class InputIterator> Kd_tree(InputIterator first, InputIterator beyond, Splitter s=Splitter(),Traits t=Traits()); 

/// @} 

/// \name Operations 
/// @{

/*! 
Inserts the point `p` in the \f$ k\f$-\f$ d\f$ tree. 
*/ 
void insert(Point_d p); 

/*! 
Inserts the elements from the sequence `[first, beyond)` in the \f$ k\f$-\f$ d\f$ tree. 
The value type of the `InputIterator` must be `Point_d`. 
*/ 
template <class InputIterator> void insert(InputIterator first, InputIterator beyond); 

/*! 
Reports the points that are approximately contained by `q`. 
The types `FuzzyQueryItem::Point_d` and `Point_d` must be equivalent. 
To use this function `Traits` must be a model of the concept `RangeSearchTraits`. 

*/ 
template <class OutputIterator, class FuzzyQueryItem> 
OutputIterator search(OutputIterator it, FuzzyQueryItem q) const; 

/*! 
Returns a const iterator to the first point in the tree. 
*/ 
iterator begin() const; 

/*! 
Returns the appropriate past-the-end const iterator. 
*/ 
iterator end() const; 

/*! 
Removes all points from the \f$ k\f$-\f$ d\f$ tree. 
*/ 
void clear(); 

/*! 
Returns the number of points that are stored in the tree. 
*/ 
size_type size() const; 

/*! 
return the instance of the traits used to construct the tree. 
*/ 
Traits traits() const; 

/*! 
  \advanced Returns a handle to the root node of the tree. 
*/ 
Node_handle root(); 

/*! 
  \advanced Returns a const handle to the root node of the tree. 
*/ 
Node_const_handle root() const; 

/*!  \advanced returns a const reference to the bounding box of the
root node of the tree.
*/ 
const Kd_tree_rectangle<FT>& bounding_box() const; 

/*! 
  \advanced Inserts statistics of the tree into the output stream `s`. 
*/ 
std::ostream& statistics(std::ostream& s) const; 

/// @}

}; /* end Kd_tree */
} /* end namespace CGAL */
