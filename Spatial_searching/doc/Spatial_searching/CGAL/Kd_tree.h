namespace CGAL {

/*!
\ingroup SearchClasses

The class `Kd_tree` defines a `k-d` tree.

\cgalHeading{Parameters}

\tparam Traits must be a model of the concept
`SearchTraits`, for example `Search_traits_2<Simple_cartesian<double> >`.

\tparam Splitter must be a model for the concept `Splitter`.
It defaults to `Sliding_midpoint<Traits>`.

\tparam UseExtendedNode must be  `Tag_true`, if the
tree shall be built with extended nodes, and `Tag_false` otherwise.

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
Dimension tag.
*/
typedef unspecified_type D;

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
typedef unspecified_type Splitter;

/*!
Bidirectional const iterator with value type `Point_d` that allows
to enumerate all points in the tree.
*/
typedef unspecified_type iterator;

/*!
A handle with value type `Kd_tree_node<Traits,Splitter>`.
*/
typedef unspecified_type Node_handle;

/*!
A const handle with value type `Kd_tree_node<Traits,Splitter>`.
*/
typedef unspecified_type Node_const_handle;

/*!
Random access const iterator with value type `const Point_d*`.
*/
typedef unspecified_type Point_d_iterator;

/*!
A type that counts the number of elements in a `k-d` tree.
*/
typedef unspecified_type size_type;

/// @}

/// \name Creation
/// @{

/*!
Constructs an empty `k-d` tree.
*/
Kd_tree(Splitter s=Splitter(),Traits t=Traits());

/*!

Constructs a `k-d` tree on the elements from the sequence
`[first, beyond)` using the splitting rule implemented by `s`.
The value type of the `InputIterator` must be `Point_d`.

*/
template <class InputIterator> Kd_tree(InputIterator first, InputIterator beyond, Splitter s=Splitter(),Traits t=Traits());

/// @}

/// \name Operations
/// @{

/*!
Inserts the point `p` in the `k-d` tree.
*/
void insert(Point_d p);

/*!
Inserts the elements from the sequence `[first, beyond)` in the `k-d` tree.
The value type of the `InputIterator` must be `Point_d`.
*/
template <class InputIterator> void insert(InputIterator first, InputIterator beyond);


/*
Pre-allocates memory in order to store at least 'size' points.
*/
void reserve(size_t size);

/*
Returns the number of points for which memory has been pre-allocated.
 */
size_t capacity();

/*!
Reports the points that are approximately contained by `q`.
The types `FuzzyQueryItem::Point_d` and `Point_d` must be equivalent.
To use this function `Traits` must be a model of the concept `RangeSearchTraits`.

*/
template <class OutputIterator, class FuzzyQueryItem>
OutputIterator search(OutputIterator it, FuzzyQueryItem q) const;

/*!
Returns a const iterator to the first point in the tree.
\note Starting with %CGAL 4.6, the order of the points in the iterator range
`[begin() , end())` is not the order of  insertion of the points into the tree.
This was not guaranteed before but might have beeen observed and exploited.
*/
iterator begin() const;

/*!
Returns the appropriate past-the-end const iterator.
*/
iterator end() const;

/*!
Removes all points from the `k-d` tree.
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
Returns a handle to the root node of the tree.
*/
Node_handle root();

/*!
Returns a const handle to the root node of the tree.
*/
Node_const_handle root() const;

/*!
Returns a const reference to the bounding box of the
root node of the tree.
*/
const Kd_tree_rectangle<FT,D>& bounding_box() const;

/*!
Inserts statistics of the tree into the output stream `s`.
*/
std::ostream& statistics(std::ostream& s) const;

/// @}

}; /* end Kd_tree */
} /* end namespace CGAL */
