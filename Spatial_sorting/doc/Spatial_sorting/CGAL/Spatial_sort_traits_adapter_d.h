namespace CGAL {

/*!
\ingroup PkgSpatialSortingTraits

Given a property map associating a key to a point, the class `Spatial_sort_traits_adapter_d` induces a spatial 
reorder of the keys instead of the points, the comparisons being done on the associated points. 
In other words, the traits provides to a spatial sort algorithm a point type which is a key, 
while the actual point type is `Base_traits::Point_d`. 


\tparam Base_traits must be a model for `SpatialSortingTraits_d`. 
\tparam PointPropertyMap must be a model of `ReadablePropertyMap`
with value type `Base_traits::Point_d`.

\cgalModels `SpatialSortingTraits_d`

*/
template< typename Base_traits, typename PointPropertyMap >
class Spatial_sort_traits_adapter_d : public Base_traits {
public:

/// \name Types 
/// @{

/*!

*/ 
boost::property_traits<PointPropertyMap>::key_type Point_d; 

/// @} 

/// \name Creation 
/// @{

/*!

*/ 
Spatial_sort_traits_adapter_d(Base_traits base=Base_traits()); 

/*!

*/ 
Spatial_sort_traits_adapter_d(const PointPropertyMap& ppmap,Base_traits base=Base_traits()); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns a const reference to the point property map. 
*/ 
const PointPropertyMap& point_property_map() const; 

/// @}

}; /* end Spatial_sort_traits_adapter_d */
} /* end namespace CGAL */
