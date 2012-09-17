namespace CGAL {

/*!
\ingroup SearchClasses

The class `Orthogonal_k_neighbor_search` implements approximate \f$ k\f$-nearest and 
\f$ k\f$-furthest neighbor searching on a tree 
using an orthogonal distance class. 

### Parameters ###

Expects for the first template argument an implementation of the concept `SearchTraits`, 
for example `CGAL::Search_traits_2<CGAL::Cartesian<double> >`. 

Expects for the second template argument a model of the 
concept `GeneralDistance`. If `Traits` is 
`CGAL::Search_traits_adapter<Key,PointPropertyMap,BaseTraits>` 
the default type is `CGAL::Distance_adapter<Key,PointPropertyMap,CGAL::Euclidean_distance<Traits> >`, 
and `CGAL::Euclidean_distance<Traits>` otherwise. 

The default type is 
`CGAL::Euclidean_distance<Traits>`. 

Expects for third template argument a model of the concept `Splitter`. 
The default type is `CGAL::Sliding_midpoint<Traits>`. 

Expects for fourth template argument an implementation of the concept `SpatialTree`. 
The default type is `CGAL::Kd_tree<Traits, Splitter, CGAL::Tag_true>`. The 
template argument must be `CGAL::Tag_true` because orthogonal search needs extended 
kd tree nodes. 

\sa `CGAL::K_neighbor_search<Traits, GeneralDistance, Splitter, SpatialTree>` 

*/
template< typename Traits, typename OrthogonalDistance, typename Splitter, typename SpatialTree >
class Orthogonal_k_neighbor_search {
public:

/// \name Types 
/// @{

/*! 
Point type. 
*/ 
Traits::Point_d Point_d; 

/*! 
Number type. 
*/ 
Traits::FT FT; 

/*! 
Distance type. 
*/ 
OrthogonalDistance Distance; 

/*! 
Query item. 
*/ 
GeneralDistance::Query_item Query_item; 

/*! 
Pair of point and transformed distance. 
*/ 
std::pair<Point_d,FT> Point_with_transformed_distance; 

/*! 
Bidirectional const iterator with value type `Point_with_transformed_distance` 
for enumerating approximate neighbors. 
*/ 
typedef Hidden_type iterator; 

/*! 
The tree type. 
*/ 
SpatialTree Tree; 

/// @} 

/// \name Operations 
/// @{

/*! 
Constructor for searching approximately \f$ k\f$ neighbors of the query item `query` 
in the points stored in `tree` using 
distance `d` and approximation factor `eps`.`sorted` indicates 
if the computed sequence of \f$ k\f$w-nearest neighbors needs to be sorted. 
*/ 
Orthogonal_k_neighbor_search(SpatialTree tree, Query_item query, unsigned int k=1, FT eps=FT(0.0), 
bool search_nearest=true, 
OrthogonalDistance d=OrthogonalDistance(),bool sorted=true); 

/*! 
Returns a const iterator to the approximate nearest or furthest neighbor. 
*/ 
iterator begin() const; 

/*! 
Returns the appropriate past-the-end const iterator. 
*/ 
iterator end() const; 

/*! 
\advanced Inserts statistics of the search process into the output stream `s`. 
*/ 
std::ostream& statistics(std::ostream& s); 

/// @}

}; /* end Orthogonal_k_neighbor_search */
} /* end namespace CGAL */
