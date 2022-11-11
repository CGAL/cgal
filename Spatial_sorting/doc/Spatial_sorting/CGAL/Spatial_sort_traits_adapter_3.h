namespace CGAL {

/*!
\ingroup PkgSpatialSortingTraits

Given a property map associating a key to a point, the class `Spatial_sort_traits_adapter_3` induces a spatial
reorder of the keys instead of the points, the comparisons being done on the associated points.
In other words, the traits provides to a spatial sort algorithm a point type which is a key,
while the actual point type is `Base_traits::Point_3`.


\tparam Base_traits must be a model for `SpatialSortingTraits_3`.
\tparam PointPropertyMap must be a model of `ReadablePropertyMap`
with value type `Base_traits::Point_3`.

\cgalModels `SpatialSortingTraits_3`

*/
template< typename Base_traits, typename PointPropertyMap >
class Spatial_sort_traits_adapter_3 : public Base_traits {
public:

/// \name Types
/// @{

/*!

*/
typedef boost::property_traits<PointPropertyMap>::key_type Point_3;

/// @}

/// \name Creation
/// @{

/*!

*/
Spatial_sort_traits_adapter_3(Base_traits base=Base_traits());

/*!

*/
Spatial_sort_traits_adapter_3(const PointPropertyMap& ppmap,Base_traits base=Base_traits());

/// @}

/// \name Operations
/// @{

/*!
Returns a const reference to the point property map.
*/
const PointPropertyMap& point_property_map() const;

/// @}

}; /* end Spatial_sort_traits_adapter_3 */
} /* end namespace CGAL */
