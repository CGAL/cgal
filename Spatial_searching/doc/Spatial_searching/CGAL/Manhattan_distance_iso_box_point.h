namespace CGAL {

/*!
\ingroup DistanceClasses

The class `Manhattan_distance_iso_box_point` provides an implementation of 
the `GeneralDistance` concept 
for the Manhattan distance (\f$ l_1\f$ metric) 
between a `d`-dimensional iso-box and a `d`-dimensional point 
and the Manhattan distance between a `d`-dimensional iso-box 
and a `d`-dimensional iso-box defined as a `k-d` tree rectangle. 

\tparam Traits must be a model for the concept 
`SearchTraits`, for example `Search_traits_3<Simple_cartesian<double> >`. 

\cgalModels `GeneralDistance`

\sa `GeneralDistance` 

*/
template< typename Traits >
class Manhattan_distance_iso_box_point {
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
typedef Traits::Iso_box_d Query_item; 

/// @} 

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
Manhattan_distance_iso_box_point(Traits t=Traits()); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the transformed 
distance between `b` and `p`. 
*/ 
FT transformed_distance(Query_item b, Point_d p) const; 

/*!
Returns the transformed value of of `d`. 
*/ 
FT transformed_distance(FT d) const; 

/*!
Returns the value of the inverse 
of the transform function applied to `d`. 
*/ 
FT inverse_of_transformed_distance(FT d) const; 

/*!
Returns the minimal distance between a point from `b` and a point from 
`r`. 
*/ 
FT min_distance_to_rectangle(Query_item b, Kd_tree_rectangle<FT,D> r) const; 

/*!
Returns the maximal distance between the iso-box `b` and 
a point from `r` furthest to `b`. 
*/ 
FT max_distance_to_rectangle(Query_item b, Kd_tree_rectangle<FT,D> r) const; 

/// @}

}; /* end Manhattan_distance_iso_box_point */
} /* end namespace CGAL */
