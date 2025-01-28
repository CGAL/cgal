namespace CGAL {

/*!
\ingroup PkgConvexHull2Traits

The class `Convex_hull_traits_adapter_2` serves as a traits class for all the two-dimensional
convex hull and extreme point calculation functions.

Given a property map associating a key to a point, the class `Convex_hull_traits_adapter_2` enables
to compute the sequence of keys for which the associated points form a convex hull,
performing the predicates of the base traits class on the points associated to the keys.

\cgalModels{ConvexHullTraits_2}

\sa `CGAL::Convex_hull_constructive_traits_2<R>`
\sa `CGAL::Projection_traits_xy_3<K>`
\sa `CGAL::Projection_traits_yz_3<K>`
\sa `CGAL::Projection_traits_xz_3<K>`

*/
template< typename BaseTraits, typename PointPropertyMap >
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
Convex_hull_traits_2(BaseTraits base=BaseTraits());

/*!
*/
Convex_hull_traits_2(const PointPropertyMap& ppmap, BaseTraits base=BaseTraits());
/// @}

/// \name Operations
/// @{

/*!
Returns a const reference to the point property map.
*/
const PointPropertyMap& point_property_map() const;


/// @}

}; /* end Convex_hull_traits_adapter_2 */
} /* end namespace CGAL */
