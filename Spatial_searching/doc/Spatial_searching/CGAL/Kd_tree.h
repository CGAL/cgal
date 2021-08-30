namespace CGAL {

/*!
\ingroup SearchClasses

The class `Kd_tree` defines a `k-d` tree.

\tparam Traits must be a model of the concept
`SearchTraits`, for example `Search_traits_2<Simple_cartesian<double> >`.

\tparam Splitter must be a model for the concept `Splitter`.
It defaults to `Sliding_midpoint<Traits>`.

\tparam UseExtendedNode must be  `Tag_true`, if the
tree shall be built with extended nodes, and `Tag_false` otherwise.

\tparam EnablePointsCache can be `Tag_true` or `Tag_false`.
Not storing the points coordinates inside the tree usually generates a
lot of cache misses, leading to non-optimal performance. This is the case
for example when indices are stored inside the tree,
or to a lesser extent when the points coordinates are stored
in a dynamically allocated array (e.g., `Epick_d` with dynamic
dimension) &mdash; we says "to a lesser extent" because the points
are re-created by the kd-tree in a cache-friendly order after its construction,
so the coordinates are more likely to be stored in a near-optimal order on the
heap. When `EnablePointsCache` is set to `Tag_true`, the points
coordinates will be cached in an optimal way. This will
increase memory consumption but provide better search performance.
See also the `GeneralDistance` and `FuzzyQueryItem` concepts for
additional requirements when using such a cache.

\sa `CGAL::Kd_tree_node<Traits>`
\sa `CGAL::Search_traits_2<Kernel>`
\sa `CGAL::Search_traits_3<Kernel>`
\sa `CGAL::Search_traits<FT_,Point,CartesianIterator,ConstructCartesianIterator>`

*/
template< typename Traits, typename Splitter, typename UseExtendedNode, typename EnablePointsCache >
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

/*!
The constructor does not build the internal data structure, and it
is also not updated after calls to `insert()`.
The method `build()` is called implicitly
at the first call to a query or removal member function. You can call
`build()` explicitly to ensure that the next call to
query functions will not trigger the reconstruction of the
data structure.

\tparam ConcurrencyTag enables sequential versus parallel
algorithm. Possible values are `Sequential_tag`, `Parallel_tag`, and
`Parallel_if_available_tag`. This template parameter is optional:
calling `build()` without specifying the concurrency tag will result
in `Sequential_tag` being used. If `build()` is not called by the user
but called implicitly at the first call to a query or removal member
function, `Sequential_tag` is also used.

*/
template <typename ConcurrencyTag>
void build();

/*!
This clears the internal data structure, which then gets rebuilt either by an
explicit call to `build()` or implicitly by the next query or removal. The only
reason to call this function explicitly is to rebalance the tree after some
number of removals.
*/
void invalidate_build();
/// @}

/// \name Operations
/// @{

/*!
Inserts the point `p` in the `k-d` tree.
\note Insertions do not dynamically update the internal data structure. The
next query, or a call to `build()`, automatically triggers a rebuild of the
whole structure.
*/
void insert(Point_d p);

/*!
Inserts the elements from the sequence `[first, beyond)` in the `k-d` tree.
The value type of the `InputIterator` must be `Point_d`.
*/
template <class InputIterator> void insert(InputIterator first, InputIterator beyond);

/*!
Removes the point `p` from the `k-d` tree. It uses `equal_to_p` to identify
the point after locating it, which can matter in particular when 2 points are
in the same place. `IdentifyPoint` is a unary functor that takes a `Point_d`
and returns a `bool`.  This is a limited and naive implementation that does not
rebalance the tree. On the other hand, the tree remains valid and ready for
queries. If the internal data structure is not already built, for instance
because the last operation was an insertion, it first calls `build()`.
*/
template<class IdentifyPoint>
void remove(Point_d p, IdentifyPoint identify_point);

/*!
Removes point `p`, calling the 2-argument function `remove()` with a functor
that simply compares coordinates.
*/
void remove(Point_d p);

/*
Pre-allocates memory in order to store at least 'size' points.
*/
void reserve(size_t size);

/*
Returns the number of points for which memory has been pre-allocated.
 */
size_t capacity();

/*!
Reports any point that is approximately contained by `q`.
The types `FuzzyQueryItem::Point_d` and `Point_d` must be equivalent.
To use this function `Traits` must be a model of the concept `RangeSearchTraits`.
*/
template <class FuzzyQueryItem>
boost::optional<Point_d> search_any_point(FuzzyQueryItem q) const;

/*!
Reports the points that are approximately contained by `q`.
The types `FuzzyQueryItem::Point_d` and `Point_d` must be equivalent.
To use this function `Traits` must be a model of the concept `RangeSearchTraits`.

*/
template <class OutputIterator, class FuzzyQueryItem>
OutputIterator search(OutputIterator it, FuzzyQueryItem q) const;

/*!
Returns a const iterator to the first point in the tree.
\note Starting with \cgal 4.6, the order of the points in the iterator range
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
