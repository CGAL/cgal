namespace CGAL {

/*!
\ingroup SearchClasses

The class `K_neighbor_search` implements approximate \f$ k\f$-nearest and \f$ k\f$-furthest neighbor searching 
using standard search on a tree using a general distance class. The tree may be built with 
extended or unextended nodes. 

Parameters 
-------------- 

Expects for the first template argument an implementation of the concept `SearchTraits`, 
for example `CGAL::Cartesian_d<double>`. 

Expects for the second template argument a model of the 
concept `GeneralDistance`. If `Traits` is 
`CGAL::Search_traits_adapter<Key,PointPropertyMap,BaseTraits>` 
the default type is `CGAL::Distance_for_point_adapter<Key,PointPropertyMap,CGAL::Euclidean_distance<Traits> >`, 
and `CGAL::Euclidean_distance<Traits>` otherwise. 

Expects for fourth template argument an implementation of the concept `SpatialTree`. 
The default type is `CGAL::Kd_tree<Traits, Splitter, CGAL::Tag_false>`. The 
template argument `CGAL::Tag_false` makes that the tree is built with unextended nodes. 

\sa `CGAL::Orthogonal_k_neighbor_search<Traits, OrthogonalDistance, Splitter, SpatialTree>` 

*/
template< typename Traits, typename GeneralDistance, typename Splitter, typename SpatialTree >
class K_neighbor_search {
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
GeneralDistance Distance; 

/*! 
Pair of point and transformed distance. 
*/ 
std::pair<Point_d,FT> Point_with_transformed_distance; 

/*! 
Bidirectional const iterator with value type `Point_with_distance` 
for enumerating approximate neighbors. 
*/ 
typedef Hidden_type iterator; 

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
Constructor for searching approximately \f$ k\f$ neighbors of the query item `q` 
in the points stored in `tree` using 
distance class \f$ d\f$ and approximation factor `eps`. `sorted` indicates 
if the computed sequence of \f$ k\f$-nearest neighbors needs to be sorted. 
*/ 
K_neighbor_search(const Tree& tree, Query_item q, unsigned int k=1, FT eps=FT(0.0), 
bool search_nearest=true, GeneralDistance d=GeneralDistance(),bool sorted=true); 

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

\advanced Inserts statistics of the search process into the output stream `s`. 

*/ 
std::ostream& statistics(std::ostream& s); 

/// @}

}; /* end K_neighbor_search */
} /* end namespace CGAL */
