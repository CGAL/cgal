namespace CGAL {

/*!
\ingroup SearchClasses

The class `K_neighbor_search` implements approximate `k`-nearest and `k`-furthest neighbor searching 
using standard search on a tree using a general distance class. The tree may be built with 
extended or unextended nodes. 


\tparam Traits must be an implementation of the concept `SearchTraits`, 
for example `Simple_cartesian<double>`. 

\tparam Splitter must be a model of the 
concept `GeneralDistance`. If `Traits` is
`Search_traits_adapter<Key,PointPropertyMap,BaseTraits>` 
the default type is `Distance_adapter<Key,PointPropertyMap,Euclidean_distance<BaseTraits> >`, 
and `Euclidean_distance<Traits>` otherwise.

\tparam SpatialTree must be an implementation of the concept `SpatialTree`. 
The default type is `Kd_tree<Traits, Splitter, Tag_false>`. The 
template argument `Tag_false` makes that the tree is built with unextended nodes. 

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
typedef Traits::Point_d Point_d; 

/*!
Number type. 
*/ 
typedef Traits::FT FT; 

/*!
Distance type. 
*/ 
typedef GeneralDistance Distance; 

/*!
Pair of point and transformed distance. 
*/ 
typedef std::pair<Point_d,FT> Point_with_transformed_distance; 

/*!
Bidirectional const iterator with value type `Point_with_distance` 
for enumerating approximate neighbors. 
*/ 
typedef unspecified_type iterator; 

/*!
Query item type. 
*/ 
typedef GeneralDistance::Query_item Query_item; 

/*!
The tree type. 
*/ 
typedef SpatialTree Tree; 

/// @} 

/// \name Creation 
/// @{

/*!
Constructor for searching approximately `k` neighbors of the query item `q` 
in the points stored in `tree` using 
distance class `d` and approximation factor `eps`. `sorted` indicates 
if the computed sequence of `k`-nearest neighbors needs to be sorted. 
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
Inserts statistics of the search process into the output stream `s`.
*/ 
std::ostream& statistics(std::ostream& s); 

/// @}

}; /* end K_neighbor_search */
} /* end namespace CGAL */
