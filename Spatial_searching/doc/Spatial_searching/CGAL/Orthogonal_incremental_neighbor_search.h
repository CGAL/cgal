namespace CGAL {

/*!
\ingroup SearchClasses

The class `Orthogonal_incremental_neighbor_search` implements incremental nearest and furthest neighbor searching on a tree. 

### Parameters ###

Expects for the first template argument a model of the concept `SearchTraits`, 
for example `CGAL::Search_traits_2<CGAL::Cartesian<double> >`. 

Expects for the second template argument a model of the 
concept `GeneralDistance`. If `Traits` is 
`CGAL::Search_traits_adapter<Key,PointPropertyMap,BaseTraits>` 
the default type is `CGAL::Distance_adapter<Key,PointPropertyMap,CGAL::Euclidean_distance<Traits> >`, 
and `CGAL::Euclidean_distance<Traits>` otherwise. 

Expects for third template argument a model of the concept `Splitter`. 
The default type is `CGAL::Sliding_midpoint<Traits>`. 

Expects for fourth template argument a model of the concept `SpatialTree`. 
The default type is `CGAL::Kd_tree<Traits, Splitter, CGAL::Tag_true>`. The 
template argument must be `CGAL::Tag_true` because orthogonal search needs extended 
kd tree nodes. 

\sa `CGAL::Incremental_neighbor_search<Traits, GeneralDistance, SpatialTree>` 

*/
template< typename Traits, typename OrthogonalDistance, typename Splitter, typename SpatialTree >
class Orthogonal_incremental_neighbor_search {
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
typedef Traits::FT FT; 

/*! 
Distance type. 
*/ 
typedef OrthogonalDistance Distance; 

/*! 
Query item. 
*/ 
typedef OrthogonalDistance::Query_item Query_item; 

/*! 
Pair of point and transformed distance. 
*/ 
typedef std::pair<Point_d,FT> Point_with_transformed_distance; 

/*! 
const input iterator with value type `Point_with_transformed_distance` 
for enumerating approximate neighbors. 
*/ 
typedef Hidden_type iterator; 

/*! 
The tree type. 
*/ 
typedef SpatialTree Tree; 

/// @} 

/// \name Creation 
/// @{

/*! 
Constructor for incremental neighbor searching of the query item `query` 
in the points stored `tree` using a distance `d` and approximation factor `eps`. 
*/ 
Orthogonal_incremental_neighbor_search(SpatialTree& tree, Query_item query, FT eps=FT(0.0), 
bool search_nearest=true, 
OrthogonalDistance d=OrthogonalDistance()); 

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

Inserts statistics of the search process into the output stream `s`. 

*/ 
std::ostream& statistics(std::ostream& s) const; 

/// @}

}; /* end Orthogonal_incremental_neighbor_search */
} /* end namespace CGAL */
