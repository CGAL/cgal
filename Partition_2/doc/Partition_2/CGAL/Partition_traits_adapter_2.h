namespace CGAL {

/*!
\ingroup PkgPartition2Ref

Given a property map associating a key to a point, the class `Partition_traits_adapter_2` induces a spatial 
reorder of the keys instead of the points, the comparisons being done on the associated points. 
In other words, the traits provides to a spatial sort algorithm a point type which is a key, 
while the actual point type is `Base_traits::Point_2`. 


\tparam Base_traits must be a model for `SpatialSortingTraits_2`. 
\tparam PointPropertyMap must be a model of `ReadablePropertyMap`
with value type `Base_traits::Point_2`.

\cgalModels `PartitionTraits_2`
*/
template< typename Base_traits, typename PointPropertyMap >
class Partition_traits_adapter_2 : public Base_traits {
public:

/// \name Types 
/// @{

/*!

*/ 
typdef boost::property_traits<PointPropertyMap>::key_type Point_2; 

/// @} 

/// \name Creation 
/// @{

/*!

*/ 
Partition_traits_adapter_2(Base_traits base=Base_traits()); 

/*!

*/ 
Partition_traits_adapter_2(const PointPropertyMap& ppmap,Base_traits base=Base_traits()); 

/// @} 

/// \name Operations 
/// @{

/*!
Returns a const reference to the point property map. 
*/ 
const PointPropertyMap& point_property_map() const; 

/// @}

}; /* end Partition_traits_adapter_2 */
} /* end namespace CGAL */
