/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalConcept

Requirements of an orthogonal distance class supporting incremental distance updates.
To optimize distance computations transformed distances are used.
E.g., for an Euclidean distance the transformed distance is the squared Euclidean distance.

\cgalRefines `GeneralDistance`

\cgalHasModel `CGAL::Euclidean_distance<Traits>`
\cgalHasModel `CGAL::Weighted_Minkowski_distance<Traits>`

*/

class OrthogonalDistance {
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

/// \name Creation
/// @{

/*!
Constructor implementing distance for `d`-dimensional points.
*/
OrthogonalDistance();

/// @}

/// \name Operations
/// @{

/*!
Returns the transformed distance between `q` and `r`.
*/
FT transformed_distance(Query_item q, Point_d r) const;

/*!
Returns the transformed distance between `q` and
the point on the boundary of `r` closest to `q`.
*/
FT min_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT,D> r) const;

/*!
Returns the transformed distance between `q` and
the point on the boundary of `r` closest to `q`.
The vector `dists` has the size of the dimension of the data points
and is filled with the distance in each dimension.
*/
  FT min_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT,D> r, std::vector<FT>& dists);

/*!
Returns the transformed distance between `q` and
the point on the boundary of `r` farthest to `q`.
*/
FT max_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT,D> r) const;

/*!
Returns the transformed distance between `q` and
the point on the boundary of `r` farthest to `q`.
The vector `dists` has the size of the dimension of the data points
and is filled with the distance in each dimension.
*/
  FT max_distance_to_rectangle(Query_item q, Kd_tree_rectangle<FT,D> r, std::vector<FT>& dists);

/*!
Returns the transformed distance.
*/
FT transformed_distance(FT d) const;

/*!
Returns the inverse of the transformed distance.
*/
FT inverse_of_transformed_distance(FT d) const;

/*!
Updates `dist` incrementally and returns the updated distance.
*/
FT new_distance(FT dist, FT old_off, FT new_off, int cutting_dimension) const;

/// @}

}; /* end OrthogonalDistance */
