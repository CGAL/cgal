namespace CGAL {

/*!
\ingroup PkgConvexHull2Traits

The class `Convex_hull_traits_adapter_2` serves as a traits class for all the two-dimensional 
convex hull and extreme point calculation function. This class corresponds 
to the default traits class for these functions. 

Given a property map associating a key to a point, the class `Convex_hull_traits_adapter_2` enables
to compute the sequence of keys for which the associted points are the convex hull points , 
the predicates being computed on the associated points. 
In other words, the traits provides to a convex hull algorithm a point type which is a key, 
while the actual point type is `Base_traits::Point_2`. 

\cgalModels `ConvexHullTraits_2`

\sa `CGAL::Convex_hull_constructive_traits_2<R>` 
\sa `CGAL::Projection_traits_xy_3<K>`
\sa `CGAL::Projection_traits_yz_3<K>`
\sa `CGAL::Projection_traits_xz_3<K>`

*/
template< typename Base_traits, typename PointPropertyMap >
class Convex_hull_traits_adapter_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef boost::property_traits<PointPropertyMap>::key_type Point_2; 

/// @} 

/// \name Creation 
/// @{

/*!
*/ 
Convex_hull_traits_2(Base_traits base=Base_traits()); 

/*!
*/ 
Convex_hull_traits_2(const PointPropertyMap& ppmap, Base_traits base=Base_traits()); 
/// @} 

/// \name Operations 
/// @{
  
/*!
Returns a const reference to the point property map. 
*/ 
const PointPropertyMap& point_property_map() const; 


/// @}

}; /* end Convex_hull_traits_2 */
} /* end namespace CGAL */
