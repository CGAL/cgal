namespace CGAL {

/*!
\ingroup PkgPartition2Ref

Given a property map associating a key to a point, the class `Partition_traits_adapter_2`
enables to compute a partition of a sequence of keys, performing the predicates of the base traits class
on the points associated to the keys.



\tparam Base_traits must be a model for `PartitionTraits_2`. 
\tparam PointPropertyMap must be a model of `ReadablePropertyMap`
with value type `Base_traits::Point_2`.

\cgalModels `PartitionTraits_2`
*/
template< typename BaseTraits, typename PointPropertyMap >
class Partition_traits_adapter_2 : public BaseTraits {
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
Partition_traits_adapter_2(BaseTraits base=BaseTraits()); 

/*!

*/ 
Partition_traits_adapter_2(const PointPropertyMap& ppmap,BaseTraits base=BaseTraits()); 

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
