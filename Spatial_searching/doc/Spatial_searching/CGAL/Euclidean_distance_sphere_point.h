namespace CGAL {

/*!
\ingroup DistanceClasses

The class `Euclidean_distance_sphere_point` provides an implementation of the 
`GeneralDistance` concept for the Euclidean distance (\f$ l_2\f$ 
metric) between a \f$ d\f$-dimensional sphere and a point, and the 
Euclidean distance between a \f$ d\f$-dimensional sphere and a 
\f$ d\f$-dimensional iso-rectangle defined as a \f$k\f$-\f$d\f$ tree rectangle. 


\tparam Traits must be a model of the concept `SearchTraits`, 
for example `Simple_cartesian_d<double>`. 

\cgalModels `GeneralDistance`

\sa `GeneralDistance` 

*/
template< typename Traits >
class Euclidean_distance_sphere_point {
public:

/// \name Types 
/// @{

/*!
Dimension Tag.
*/
typedef unspecified_type D;

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
typedef Traits::Sphere_d Sphere_d; 

/// @} 

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Euclidean_distance_sphere_point(Traits t=Traits()); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the distance between `s` and `p`. 
*/ 
NT transformed_distance(Query_item s, Point_d p) const; 

/*!
Returns the minimal distance between a point from the sphere `s` and a point from 
`r`. 
*/ 
NT min_distance_to_rectangle(Query_item s, Kd_tree_rectangle<FT,D> r) const; 

/*!
Returns the maximal distance between the sphere `s` and 
a point from `r` furthest to `s`. 
*/ 
NT max_distance_to_rectangle(Query_item s, Kd_tree_rectangle<FT,D> r) const; 

/*!
Returns \f$ d^2\f$. 
*/ 
NT transformed_distance(NT d) const; 

/*!
Returns \f$ d^{1/2}\f$. 
*/ 
NT inverse_of_transformed_distance(NT d) const; 

/// @}

}; /* end Euclidean_distance_sphere_point */
} /* end namespace CGAL */
