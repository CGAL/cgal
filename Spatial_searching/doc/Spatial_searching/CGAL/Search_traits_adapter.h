namespace CGAL {

/*!
\ingroup DistanceClasses

A class that uses a point property map to adapt a distance class to work on a
key as point type.
When using `Search_traits_adapter<Key,PointPropertyMap,BaseTraits>`
in a nearest neighbor search algorithm, this class must be used as distance.


\tparam Key is a type that is associated to a point of type `Base_distance::Point_d`.

\tparam PointPropertyMap is a model of `ReadablePropertyMap`
with `Key` as  key type  and `Base_distance::Point_d` as  value type.

\tparam Base_distance is a model of either `GeneralDistance` or `OrthogonalDistance`.

\cgalModels{GeneralDistance if `Base_distance` is a model of `GeneralDistance`,
            OrthogonalDistance if `Base_distance` is a model of `OrthogonalDistance`}

\sa `Search_traits_adapter<Key,PointPropertyMap,BaseTraits>`

*/
template< typename Key, typename PointPropertyMap, typename Base_distance >
class Distance_adapter : Base_distance {
public:

/// \name Types
/// @{

/*!

*/
typedef Base_distance::FT FT;

/*!

*/
typedef Key Point_d;

/*!

*/
typedef Base_distance::Query_item Query_item;

/// @}

/// \name Creation
/// @{

/*!
Constructor initializing the class to `base` and setting the point property map of the class to `ppmap`.
*/
Distance_adapter(const PointPropertyMap& ppmap=PointPropertyMap(),const Base_distance& base=Base_distance());

/// @}

/// \name Operations
/// @{

/*!
Returns the point property map.
*/
const PointPropertyMap& point_property_map() const;

/// @}

}; /* end Distance_adapter */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup SearchTraitsClasses

The class `Search_traits_adapter` can be used as a template parameter of the kd tree
and the search classes. When using this traits class in a nearest neighbor search algorithm, the class
`Distance_adapter<Key,PointPropertyMap,Base_distance>`
must be used as distance.


\tparam Key is a type that is associated to a point of type `Base_distance::Point_d`.

\tparam PointPropertyMap is a model of `LvaluePropertyMap`
with `Key` as `key_type` and `Base_distance::Point_d` as `value_type`.

\tparam BaseTraits is a model of either `SearchTraits` or `RangeSearchTraits`.

\cgalModels{SearchTraits if `BaseTraits` is a model of `SearchTraits`.,
            RangeSearchTraits if `BaseTraits` is a model of `RangeSearchTraits`.}

\sa `Distance_adapter<Key,PointPropertyMap,Base_distance>`
\sa `Search_traits_2<Kernel>`
\sa `Search_traits_3<Kernel>`
\sa `Search_traits_d<Kernel>`
\sa `Search_traits<Point,CartesianConstIterator,ConstructCartesianConstIterator>`

*/
template< typename Key, typename PointPropertyMap, typename BaseTraits >
class Search_traits_adapter : public BaseTraits {
public:

/// \name Types
/// @{

/*!

*/
typedef BaseTraits::Dimension Dimension;

/*!

*/
typedef Key Point_d;

/*!

*/
typedef BaseTraits::FT FT;

/*!

*/
typedef BaseTraits::Cartesian_const_iterator_d Cartesian_const_iterator_d;

/*!

*/
typedef BaseTraits Base;

/// @}

/// \name Creation
/// @{

/*!
Constructor initializing the class to `base` and setting the point property map of the class to `ppmap`.
*/
Search_traits_adapter(const PointPropertyMap& ppmap=PointPropertyMap(),const BaseTraits& base=BaseTraits());

/// @}

/// \name Operations
/// @{

/*!
Returns the point property map.
*/
const PointPropertyMap& point_property_map() const;

/// @}

}; /* end Search_traits_adapter */
} /* end namespace CGAL */
