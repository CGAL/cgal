namespace CGAL {

/*!
\ingroup DistanceClasses

The class `Euclidean_distance` provides an implementation of the concept `OrthogonalDistance`, with the 
Euclidean distance (\f$ l_2\f$ metric). 
To optimize distance computations squared distances are used. 

\tparam Traits must be a model of the concept 
`SearchTraits`, for example `Search_traits_2<Simple_cartesian<double> >`. 

\cgalModels `OrthogonalDistance`

\sa `OrthogonalDistance` 
\sa `CGAL::Weighted_Minkowski_distance<Traits>` 

*/
template< typename Traits >
class Euclidean_distance {
public:

/// \name Types 
/// @{

/*!
Dimension Tag.
*/
typedef Traits::Dimension D;


/*!
Number type. 
*/ 
typedef Traits::FT FT; 

/*!
Point type. 
*/ 
typedef Traits::Point_d Point_d; 

/*!
Query item type. 
*/ 
typedef Point_d Query_item; 

/// @} 

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Euclidean_distance(Traits t=Traits()); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the squared Euclidean distance between `q` and `p`. 
*/ 
FT transformed_distance(Query_item q, Point_d p) const; 

/*!
Returns the squared Euclidean distance between `q` and 
the point on the boundary of `r` closest to `q`. 
*/ 
FT min_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT,D> r) const; 

/*!
Returns the squared Euclidean distance between `q` and 
the point on the boundary of `r` closest to `q`. Stores 
the distances in each dimension in `dists`. 
*/ 
FT min_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT,D> r, vector<FT>& dists); 

/*!
Returns the squared Euclidean distance, where \f$ d\f$ denotes the distance between `q` and 
the point on the boundary of `r` farthest to `q`. 
*/ 
FT max_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT,D> r) const; 

/*!
Returns the squared Euclidean distance, where \f$ d\f$ denotes the distance between `q` and 
the point on the boundary of `r` farthest to `q`. Stores the distances in 
each dimension in `dists`.
*/ 
FT max_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT,D> r, vector<FT>& dists); 

/*!
Updates the squared `dist` incrementally 
and returns the updated squared distance. 
*/ 
FT new_distance(FT dist, FT old_off, FT new_off, int cutting_dimension) const; 

/*!
Returns \f$ d^2\f$. 
*/ 
FT transformed_distance(FT d) const; 

/*!
Returns \f$ d^{1/2}\f$. 
*/ 
FT inverse_of_transformed_distance(FT d) const; 

/// @}

}; /* end Euclidean_distance */
} /* end namespace CGAL */
