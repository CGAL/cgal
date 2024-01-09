/*!
\ingroup PkgMesh3Concepts
\cgalConcept

The concept `MeshDomainWithFeatures_3` refines the concept `MeshDomain_3`.
While the concept
`MeshDomain_3` only exposes the 2-dimensional and 3-dimensional features of
the domain through different queries, the concept `MeshDomainWithFeatures_3` also exposes 0 and
1-dimensional features. The exposed features of the domain are respectively called
subdomains, surface patches, curves
and corners according to their respective dimensions 3, 2, 1, and 0.

Each curve is assumed to be bounded, with only one connected component, and
without auto-intersections. Each curve is also assumed to be
oriented. Therefore it is possible to define the signed geodesic distance
between two ordered points on the same curve.

\cgalRefines{MeshDomain_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Mesh_domain_with_polyline_features_3<MD>}
\cgalHasModels{CGAL::Polyhedral_mesh_domain_with_features_3<IGT>}
\cgalHasModelsEnd

\sa `MeshDomain_3`

*/

class MeshDomainWithFeatures_3 {
public:

/// \name Types
/// @{

/*!
A type to distinguish
`MeshDomain_3` models from `MeshDomainWithFeatures_3` models.
*/
typedef CGAL::Tag_true Has_features;

/*!
Numerical type.
*/
typedef unspecified_type FT;

/*!
Point type.
*/
typedef unspecified_type Point_3;

/*!
Type of indices for curves (i.e., \f$ 1\f$-dimensional features)
of the input domain.
Must be a model of CopyConstructible, Assignable, DefaultConstructible and
LessThanComparable. The default constructed value must be the value of an edge which
does not approximate a 1-dimensional feature of the input domain.
*/
typedef unspecified_type Curve_index;

/*!
Type of indices for corners (i.e., \f$ 0\f$--dimensional features)
of the input domain.
Must be a model of CopyConstructible, Assignable, DefaultConstructible and
LessThanComparable.
*/
typedef unspecified_type Corner_index;

/// @}

/*! \name Operations

*/
/// @{

/*!

Returns a point on the curve with index `ci`
at signed geodesic distance `d` from point `p`.
\pre Point `p` is supposed to be on curve `ci`. If `d > 0`, the signed
geodesic distance from `p` to the endpoint of `ci` should be greater than
`d`. If ` d < 0`, the signed geodesic distance from `p` to the origin
of the curve should be less than `d`.

*/
Point_3 construct_point_on_curve(
const Point_3& p, const Curve_index& ci, FT d) const;

/// @}

/// \name Queries
/// @{

/*!
Returns the length of the curve segment from `p` to `q`, on the curve
with index `curve_index`.

If the curve with index `curve_index` is a loop, the
orientation identifies which portion of the loop corresponds to the curve
segment, otherwise `orientation` must be compatible with the orientation
of `p` and `q` on the curve.
*/
FT curve_segment_length(const Point_3& p, const Point_3& q,
                        const Curve_index& curve_index,
                        CGAL::Orientation orientation) const;
/*!

Returns `CGAL::POSITIVE` if the signed geodesic distance from
`p` to `q` on the way through `r` along loop with index `ci`
is positive, `CGAL::NEGATIVE` if the distance is negative.
\pre `p != q && p != r && r != q`
*/
CGAL::Sign distance_sign_along_loop(const Point_3& p, const Point_3& q,
const Point_3& r, const Curve_index& ci) const;

/*!
Returns the sign of the geodesic distance from `p` to `q`, on the curve
with index `ci`.
If the curve with index `ci` is a loop, the function `distance_sign_along_loop()`
must be used instead.
*/
CGAL::Sign distance_sign(const Point_3& p, const Point_3& q,
                         const Curve_index& ci) const;

/*!
Returns the length of curve with index
`curve_index`
*/
FT curve_length(const Curve_index& curve_index) const;
/*!
Returns `true` if the portion of the curve of index `index`,
between the points `c1` and `c2`, is covered by the spheres of
centers `c1` and `c2` and squared radii `sq_r1` and `sq_r2`
respectively. The points `c1` and `c2` are assumed to lie on the curve.
*/
bool is_curve_segment_covered(const Curve_index& index,
                              CGAL::Orientation orientation,
                              const Point_3& c1, const Point_3& c2,
                              const FT sq_r1, const FT sq_r2) const;
/*!

Returns `true` if the curve
`ci` is a loop.
*/
bool is_loop(const Curve_index& ci) const;

/// @}

/// \name Retrieval of the input features
/// @{

/*!
Fills `corners` with the corners of the input domain.
The value type of `corners` must be `std::pair<Corner_index,Point_3>`.
*/
template <typename OutputIterator>
OutputIterator
get_corners(OutputIterator corners) const;

/*!

Fills `curves` with the curves
of the input domain.
`curves` value type must be
`std::tuple<Curve_index,std::pair<Point_3,%Index>,std::pair<Point_3,%Index> >`.
If the curve corresponding to an entry
in `curves` is not a loop, the pair of associated points should
belong to
two corners incident to the curve.
If it is a loop, then the same `Point_3` should be given twice and must be any
point on the loop.
The `%Index` values associated to the points are their indices w.r.t.\ their dimension.
*/
template <typename OutputIterator>
OutputIterator
get_curves(OutputIterator curves) const;

/// @}

/// \name Indices converters
/// @{

/*!

Returns the index to be stored at a vertex lying on the curve identified
by `curve_index`.
*/
Index index_from_curve_index(const Curve_index& curve_index) const;

/*!

Returns the `Curve_index` of the curve where lies a vertex with
dimension 1 and index `index`.
*/
Curve_index curve_index(const Index& index) const;

/*!

Returns the index to be stored at a vertex lying on the corner identified
by `corner_index`.
*/
Index index_from_corner_index(const Corner_index& corner_index) const;

/*!

Returns the `Corner_index` of the corner where lies a vertex with
dimension 0 and index `index`.
*/
Corner_index corner_index(const Index& index) const;

/// @}

}; /* end MeshDomainWithFeatures_3 */
