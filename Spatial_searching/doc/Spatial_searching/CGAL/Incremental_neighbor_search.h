namespace CGAL {

/*!
\ingroup SearchClasses

The class `Incremental_neighbor_search` implements incremental nearest and furthest neighbor searching
on a tree. The tree may have extended or unextended nodes.


\tparam Traits must be a model of the concept `SearchTraits`,
for example `Search_traits_2<Simple_cartesian<double> >`.

\tparam GeneralDistance must be a model of the
concept `GeneralDistance`. If `Traits` is
`Search_traits_adapter<Key,PointPropertyMap,BaseTraits>`
the default type is `Distance_adapter<Key,PointPropertyMap,Euclidean_distance<BaseTraits> >`,
and `Euclidean_distance<Traits>` otherwise.

\tparam Splitter must be a model of the concept `Splitter`.
The default type is `Sliding_midpoint<Traits>`.

\tparam SpatialTree must be a model of the concept `SpatialTree`.
The default type is `Kd_tree<Traits, Splitter, Tag_false, Tag_false>`. The
third template argument `Tag_false` makes that the tree is built with unextended nodes.

\sa `CGAL::Orthogonal_incremental_neighbor_search<Traits, OrthogonalDistance, Splitter, SpatialTree>`

*/
template< typename Traits, typename GeneralDistance, typename Splitter, typename SpatialTree >
class Incremental_neighbor_search {
public:

/// \name Types
/// @{

/*!
Point type.
*/
typedef Traits::Point_d Point_d;

/*!
Number type.
*/
typedef Traits::NT NT;

/*!
Distance type.
*/
typedef GeneralDistance Distance;

/*!
Pair of point and transformed distance.
*/
typedef std::pair<Point_d,NT> Point_with_transformed_distance;

/*!
const input iterator with value type `Point_with_transformed_distance`
for enumerating approximate neighbors.
*/
typedef unspecified_type iterator;

/*!
Query item type.
*/
GeneralDistance::Query_item Query_item;

/*!
The tree type.
*/
SpatialTree Tree;

/// @}

/// \name Creation
/// @{

/*!
Constructor for incremental neighbor searching of the query item `q`
in the points stored `tree` using a distance `d` and approximation factor `eps`.
*/
Incremental_neighbor_search(Tree& tree, QueryItem q, NT eps=NT(0.0),
bool search_nearest=true, GeneralDistance d=GeneralDistance());

/// @}

/// \name Operations
/// @{

/*!
Returns a const iterator to the approximate nearest or furthest neighbor.
*/
iterator begin() const;

/*!
Returns the appropriate past-the-end const iterator.
*/
iterator end() const;

/*!
Inserts statistics of the search process into the output
stream `s`.

*/
std::ostream& statistics(std::ostream& s);

/// @}

}; /* end Incremental_neighbor_search */
} /* end namespace CGAL */
