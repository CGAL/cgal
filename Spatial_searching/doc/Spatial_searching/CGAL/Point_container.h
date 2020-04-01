namespace CGAL {

/*!
\ingroup AdvancedClasses

A custom container for points used to build a tree.

Each point container holds the points from a rectangle associated with a
node of the tree.  In the remainder of this reference page this
rectangle is called the associated rectangle.  Provides a method to
split a container and a number of methods to support the
implementation of splitting rules.

\tparam Traits must be model of the concept `SearchTraits`.


\sa `SearchTraits`
\sa `SpatialSeparator`

*/
template< typename Traits >
class Point_container {
public:

/// \name Types
/// @{

/*!
Number type.
*/
typedef Traits::FT FT;

/*!
Point type.
*/
typedef Traits::Point_d Point_d;

/*!
An iterator with value type `Point_d*`.
*/
typedef unspecified_type iterator;

/*!
A const iterator with value type `const Point_d*`.
*/
typedef unspecified_type const_iterator;

/// @}

/// \name Creation
/// @{

/*!

Construct an empty container for storing `d`-dimensional points.

*/
Point_container(int d);

/*!

Construct the container of `d`-dimensional points of type `Point_d`
given by the iterator sequence `[begin, end)`.
*/
template <class InputIterator>
Point_container(int d, InputIterator begin, InputIterator end);

/// @}

/// \name Operations
/// @{

/*!
Given an empty container `c2` with the same dimension as `c`, splits `c` into
`c`and `c2` using the separator `sep`. If sliding is `true` after splitting
each container contains at least one point. Container `c` should contain at least two points.
*/
template <class SpatialSeparator>
void split(Point_container<Traits> &c2, SpatialSeparator sep, bool sliding=false);

/*!
Swap the contents of `c` and `c2`
*/
void swap(Point_container<Traits> &c2);

/*!
Recompute the bounding box of the points in the container.
*/
void recompute_tight_bounding_box();

/*!

Returns an iterator to a pointer to the first point.

*/
iterator begin();

/*!

Returns the appropriate past-the-end iterator.

*/
iterator end();

/*!

Returns a const iterator to a pointer to the first point.

*/
const_iterator begin() const;

/*!

Returns the appropriate past-the-end const iterator.

*/
const_iterator end() const;

/*!

Returns the dimension.

*/
int dimension() const;

/*!

Returns coordinate for which the pointer list is built.

*/
int built_coordinate() const;

/*!

Returns coordinate where the associated rectangle has maximal span.

*/
int max_span_coord() const;

/*!

Returns coordinate where the point coordinates have maximal span.

*/
int max_tight_span_coord() const;

/*!

Returns lower value of the interval corresponding to
`max_span_coord()`.

*/
FT max_span_lower() const;

/*!

Returns lower value of the interval corresponding to
`max_tight_span_coord()`. That is, the smallest
`max_tight_span_coord()`-th coordinate of the points in
`c`.

*/
FT max_tight_span_lower() const;

/*!

Returns upper value of the interval corresponding to
`max_span_coord()`.

*/
FT max_span_upper() const;

/*!

Returns upper value of the interval over all dimensions
without taking dimension `d` into account.

*/
FT max_span_upper_without_dim(int d) const;

/*!

Returns upper value of the interval corresponding to
`max_tight_span_coord()`.

*/
FT max_tight_span_upper() const;

/*!

Returns the size of the interval corresponding to `max_span_coord()`.

*/
FT max_spread() const;

/*!

Returns the size of the interval corresponding to `max_tight_span_coord()`.

*/
FT max_tight_spread() const;

/*!

Returns the median value of the points stored in the container for
dimension `split_coord`.

*/
FT median(int split_coord) const;

/*!
Returns the associated rectangle.
*/
const Kd_tree_rectangle<Traits> & bounding_box() const;

/*!
Returns the bounding box of the items in associated rectangle.
*/
const Kd_tree_rectangle<Traits> & tight_bounding_box();

/*!
Returns the dimension with the maximal point spread, for which after fair splitting
the ratio of the length of the longest side and the smallest side of the bounding box of
the items in associated rectangle,
does not exceed `aspect_ratio`.
*/
int max_tight_span_coord_balanced(FT aspect_ratio) const;

/*!
Returns the splitting value for fair splitting.
*/
FT balanced_fair(int d, FT aspect_ratio);

/*!
Returns the splitting value for sliding fair splitting.
*/
FT balanced_sliding_fair(int d, FT aspect_ratio);

/*!

Returns the number of points stored.

*/
std::size_t size() const;

/*!

Returns true if no points are present, false otherwise.

*/
bool empty() const;

/// @}

}; /* end Point_container */

/*!
Prints the point container `c` to the output stream `s` and returns `s`.
\relates Point_container
*/
template<class Traits>
std::ostream& operator<<(std::ostream& s, Point_container<Traits> c);

} /* end namespace CGAL */
