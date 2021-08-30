/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalConcept

Requirements of a distance class defining a distance between a query item
denoting a spatial object and a point.
To optimize distance computations transformed distances are used,
e.g., for a Euclidean distance the transformed distance is the squared
Euclidean distance.

\cgalHasModel `CGAL::Manhattan_distance_iso_box_point<Traits>`
\cgalHasModel `CGAL::Euclidean_distance_sphere_point<Traits>`

*/

class GeneralDistance {
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
typedef unspecified_type FT;

/*!
Point type.
*/
typedef unspecified_type Point_d;

/*!
Query item type.
*/
typedef unspecified_type Query_item;

/// @}

/// \name Operations
/// @{

/*!
Returns the transformed distance between `q` and `r`.
*/
FT transformed_distance(Query_item q, Point_d r);

/*!
Optional: must be defined when used with a `Kd_tree` where `EnablePointsCache`
is set to `Tag_true`.

Returns the transformed distance between `q` and the point whose Cartesian
coordinates are contained in the range [`begin`, `end`).
*/
template <typename Coord_iterator>
FT transformed_distance_from_coordinates(
  Query_item q, Coord_iterator begin, Coord_iterator end) const;

/*!
Optional: in most cases (e.g., Euclidean distance), the distance computation
algorithm knows before its end that the distance will be greater than or equal
to some given value. In this function, the computation can be stopped when the
distance is going to be greater than or equal to `stop_if_geq_to_this`. In this case,
the only requirement of the return value it to be \f$ \geq \f$ `stop_if_geq_to_this`.
Note that points cache does not have to be activated to enable this optimization.

Returns the transformed distance between `q` and the point whose Cartesian
coordinates are contained in the range [`begin`, `end`), or any value
\f$ \geq \f$ `stop_if_geq_to_this` if the transformed distance is
\f$ \geq \f$ `stop_if_geq_to_this`.
*/
template <typename Coord_iterator>
FT interruptible_transformed_distance(
  Query_item q, Coord_iterator begin, Coord_iterator end, FT stop_if_geq_to_this) const;

/*!
Returns the transformed distance between `q` and
the point on the boundary of `r` closest to `q`.
*/
FT min_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT,D> r) const;

/*!
Returns the transformed distance between `q` and
the point on the boundary of `r` furthest to `q`.
*/
FT max_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT,D> r) const;

/*!
Returns the transformed distance.
*/
FT transformed_distance(FT d) const;

/*!
Returns the inverse of the transformed distance.
*/
FT inverse_of_transformed_distance(FT d) const;

/// @}

}; /* end GeneralDistance */
