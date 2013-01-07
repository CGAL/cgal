namespace CGAL {

/*!
\ingroup SearchClasses

The class `Orthogonal_k_neighbor_search` implements approximate`k`-nearest and 
`k`-furthest neighbor searching on a tree 
using an orthogonal distance class. 

\cgalHeading{Parameters}

Expects for the first template argument an implementation of the concept `SearchTraits`, 
for example `Search_traits_2<Simple_cartesian<double> >`. 

Expects for the second template argument a model of the 
concept `GeneralDistance`. If `Traits` is 
`Search_traits_adapter<Key,PointPropertyMap,BaseTraits>` 
the default type is `Distance_adapter<Key,PointPropertyMap,Euclidean_distance<Traits> >`, 
and `Euclidean_distance<Traits>` otherwise. 

The default type is `Euclidean_distance<Traits>`. 

Expects for third template argument a model of the concept `Splitter`. 
The default type is `Sliding_midpoint<Traits>`. 

Expects for fourth template argument an implementation of the concept `SpatialTree`. 
The default type is `Kd_tree<Traits, Splitter, Tag_true>`. The 
template argument must be `Tag_true` because orthogonal search needs extended 
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
typedef GeneralDistance::Query_item Query_item; 

/*! 
Pair of point and transformed distance. 
*/ 
typedef std::pair<Point_d,FT> Point_with_transformed_distance; 

/*! 
Bidirectional const iterator with value type `Point_with_transformed_distance` 
for enumerating approximate neighbors. 
*/ 
typedef Hidden_type iterator; 

/*! 
The tree type. 
*/ 
typedef SpatialTree Tree; 

/// @} 

/// \name Operations 
/// @{

/*! 
Constructor for searching approximately `k` neighbors of the query item `query` 
in the points stored in `tree` using 
distance `d` and approximation factor `eps`. `sorted` indicates 
if the computed sequence of `k`-nearest neighbors needs to be sorted. 
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
\cgalAdvanced Inserts statistics of the search process into the output stream `s`. 
*/ 
std::ostream& statistics(std::ostream& s); 

/// @}

}; /* end Orthogonal_k_neighbor_search */
} /* end namespace CGAL */
