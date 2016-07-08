namespace CGAL {

/*!
\defgroup angle_grp CGAL::angle()
\ingroup kernel_global_function
*/
/// @{

/*!
returns `CGAL::OBTUSE`, `CGAL::RIGHT` or `CGAL::ACUTE` depending
on the angle formed by the two vectors `u` and `v`.
*/
template <typename Kernel>
Angle angle(const CGAL::Vector_2<Kernel>&u, 
const CGAL::Vector_2<Kernel>&v);

/*!

returns `CGAL::OBTUSE`, `CGAL::RIGHT` or `CGAL::ACUTE` depending
on the angle formed by the three points `p`, `q`, `r` (`q` being the vertex of
the angle). The returned value is the same as `angle(p - q, r - q)`.
*/
template <typename Kernel>
Angle angle(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r);

/*!

returns `CGAL::OBTUSE`, `CGAL::RIGHT` or `CGAL::ACUTE` depending
on the angle formed by the two vectors `pq`, `rs`. The returned value is
the same as `angle(q - p, s - r)`.
*/
template <typename Kernel>
Angle angle(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r,
const CGAL::Point_2<Kernel>& s);

/*!

returns `CGAL::OBTUSE`, `CGAL::RIGHT` or `CGAL::ACUTE` depending
on the angle formed by the three points `p`, `q`, `r` (`q` being the vertex of
the angle).
*/
template <typename Kernel>
Angle angle(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r);

/// @}

/// \defgroup area_grp CGAL::area()
/// \ingroup kernel_global_function
/// @{

/*!
returns the signed area of the triangle defined by the points `p`,
`q` and `r`. 
*/
template <typename Kernel>
Kernel::FT area(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r);

/// @}

/// \defgroup are_ordered_along_line_grp CGAL::are_ordered_along_line()
/// \ingroup kernel_global_function
/// \sa `are_strictly_ordered_along_line_grp`
/// \sa `collinear_are_ordered_along_line_grp`
/// \sa `collinear_are_strictly_ordered_along_line_grp`
/// @{

/*!

returns `true`, iff the three points are collinear and 
`q` lies between `p` and `r`.
Note that `true` is returned, if `q==p` or
`q==r`.
*/
template <typename Kernel>
bool are_ordered_along_line(const CGAL::Point_2<Kernel> &p, 
const CGAL::Point_2<Kernel> &q, 
const CGAL::Point_2<Kernel> &r);

/*!
returns `true`, iff the three points are collinear and 
`q` lies between `p` and `r`.
Note that `true` is returned, if `q==p` or
`q==r`.
*/
template <typename Kernel>
bool are_ordered_along_line(const CGAL::Point_3<Kernel> &p, 
const CGAL::Point_3<Kernel> &q, 
const CGAL::Point_3<Kernel> &r);


/// @}

/// \defgroup are_strictly_ordered_along_line_grp CGAL::are_strictly_ordered_along_line()
/// \ingroup kernel_global_function
/// \sa `are_ordered_along_line_grp`
/// \sa `collinear_are_ordered_along_line_grp`
/// \sa `collinear_are_strictly_ordered_along_line_grp`
/// @{

/*!
returns `true`, iff the three points are collinear and 
`q` lies strictly between `p` and `r`.
Note that `false` is returned, if `q==p` or
`q==r`.
*/
template <typename Kernel>
bool are_strictly_ordered_along_line(const CGAL::Point_2<Kernel> &p, 
const CGAL::Point_2<Kernel> &q, 
const CGAL::Point_2<Kernel> &r);

/*!
returns `true`, iff the three points are collinear and 
`q` lies strictly between `p` and `r`.
Note that `false` is returned, if `q==p` or
`q==r`.
*/
template <typename Kernel>
bool are_strictly_ordered_along_line(const CGAL::Point_3<Kernel> &p, 
const CGAL::Point_3<Kernel> &q, 
const CGAL::Point_3<Kernel> &r);

/// @}

/// \defgroup barycenter_grp CGAL::barycenter()
/// \ingroup kernel_global_function
/// \sa \link centroid `CGAL::centroid()`  \endlink
/// @{

/*!
compute the barycenter of the points `p1` and `p2` with corresponding
weights `w1` and `1-w1`.
*/
template <typename Kernel>
CGAL::Point_2<Kernel>
barycenter( const CGAL::Point_2<Kernel>& p1, const Kernel::FT&w1,
const CGAL::Point_2<Kernel>& p2);

/*!
compute the barycenter of the points `p1` and `p2` with corresponding
weights `w1` and `w2`. \pre `w1+w2 != 0`.
*/
template <typename Kernel>
CGAL::Point_2<Kernel>
barycenter( const CGAL::Point_2<Kernel>& p1, const Kernel::FT&w1,
const CGAL::Point_2<Kernel>& p2, const Kernel::FT&w2);

/*!
compute the barycenter of the points `p1`, `p2` and `p3` with corresponding
weights `w1`, `w2` and `1-w1-w2`.
*/
template <typename Kernel>
CGAL::Point_2<Kernel>
barycenter( const CGAL::Point_2<Kernel>& p1, const Kernel::FT&w1,
const CGAL::Point_2<Kernel>& p2, const Kernel::FT&w2,
const CGAL::Point_2<Kernel>& p3);

/*!
compute the barycenter of the points `p1`, `p2` and `p3` with corresponding
weights `w1`, `w2` and `w3`. \pre `w1+w2+w3 != 0`.
*/
template <typename Kernel>
CGAL::Point_2<Kernel>
barycenter( const CGAL::Point_2<Kernel>& p1, const Kernel::FT&w1,
const CGAL::Point_2<Kernel>& p2, const Kernel::FT&w2,
const CGAL::Point_2<Kernel>& p3, const Kernel::FT&w3);

/*!
compute the barycenter of the points `p1`, `p2`, `p3` and `p4` with corresponding
weights `w1`, `w2`, `w3` and `1-w1-w2-w3`.
*/
template <typename Kernel>
CGAL::Point_2<Kernel>
barycenter( const CGAL::Point_2<Kernel>& p1, const Kernel::FT&w1,
const CGAL::Point_2<Kernel>& p2, const Kernel::FT&w2,
const CGAL::Point_2<Kernel>& p3, const Kernel::FT&w3,
const CGAL::Point_2<Kernel>& p4);

/*!
compute the barycenter of the points `p1`, `p2`, `p3` and `p4` with corresponding
weights `w1`, `w2`, `w3` and `w4`. \pre `w1+w2+w3+w4 != 0`.
*/
template <typename Kernel>
CGAL::Point_2<Kernel>
barycenter( const CGAL::Point_2<Kernel>& p1, const Kernel::FT&w1,
const CGAL::Point_2<Kernel>& p2, const Kernel::FT&w2,
const CGAL::Point_2<Kernel>& p3, const Kernel::FT&w3,
const CGAL::Point_2<Kernel>& p4, const Kernel::FT&w4);

/*!
compute the barycenter of the points `p1` and `p2` with corresponding
weights `w1` and `1-w1`.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
barycenter( const CGAL::Point_3<Kernel>& p1, const Kernel::FT&w1,
const CGAL::Point_3<Kernel>& p2);

/*!
compute the barycenter of the points `p1` and `p2` with corresponding
weights `w1` and `w2`. \pre `w1+w2 != 0`.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
barycenter( const CGAL::Point_3<Kernel>& p1, const Kernel::FT&w1,
const CGAL::Point_3<Kernel>& p2, const Kernel::FT&w2);

/*!
compute the barycenter of the points `p1`, `p2` and `p3` with corresponding
weights `w1`, `w2` and `1-w1-w2`.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
barycenter( const CGAL::Point_3<Kernel>& p1, const Kernel::FT&w1,
const CGAL::Point_3<Kernel>& p2, const Kernel::FT&w2,
const CGAL::Point_3<Kernel>& p3);

/*!
compute the barycenter of the points `p1`, `p2` and `p3` with corresponding
weights `w1`, `w2` and `w3`. \pre `w1+w2+w3 != 0`.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
barycenter( const CGAL::Point_3<Kernel>& p1, const Kernel::FT&w1,
const CGAL::Point_3<Kernel>& p2, const Kernel::FT&w2,
const CGAL::Point_3<Kernel>& p3, const Kernel::FT&w3);

/*!
compute the barycenter of the points `p1`, `p2`, `p3` and `p4` with corresponding
weights `w1`, `w2`, `w3` and `1-w1-w2-w3`.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
barycenter( const CGAL::Point_3<Kernel>& p1, const Kernel::FT&w1,
const CGAL::Point_3<Kernel>& p2, const Kernel::FT&w2,
const CGAL::Point_3<Kernel>& p3, const Kernel::FT&w3,
const CGAL::Point_3<Kernel>& p4);

/*!
compute the barycenter of the points `p1`, `p2`, `p3` and `p4` with corresponding
weights `w1`, `w2`, `w3` and `w4`. \pre `w1+w2+w3+w4 != 0`.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
barycenter( const CGAL::Point_3<Kernel>& p1, const Kernel::FT&w1,
const CGAL::Point_3<Kernel>& p2, const Kernel::FT&w2,
const CGAL::Point_3<Kernel>& p3, const Kernel::FT&w3,
const CGAL::Point_3<Kernel>& p4, const Kernel::FT&w4);

/// @}

/// \defgroup bisector_grp CGAL::bisector()
/// \ingroup kernel_global_function
/// @{

/*!
constructs the bisector line of the two points `p` and `q`.
The bisector is oriented in such a way that `p` lies on its
positive side. \pre `p` and `q` are not equal.
*/
template <typename Kernel>
CGAL::Line_2<Kernel> bisector(const CGAL::Point_2<Kernel> &p,
const CGAL::Point_2<Kernel> &q);

/*!
constructs the bisector of the two lines `l1` and `l2`.
In the general case, the bisector has the direction of the vector which
is the sum of the normalized directions of the two lines, and which passes
through the intersection of `l1` and `l2`.
If `l1` and `l2` are parallel, then the bisector is defined as the line
which has the same direction as `l1`, and which is at the same distance
from `l1` and `l2`.
This function requires that `Kernel::RT` supports the `sqrt()`
operation.
*/
template <typename Kernel>
CGAL::Line_2<Kernel> bisector(const CGAL::Line_2<Kernel> &l1,
const CGAL::Line_2<Kernel> &l2);

/*!
constructs the bisector plane of the two points `p` and `q`.
The bisector is oriented in such a way that `p` lies on its
positive side. \pre `p` and `q` are not equal.
*/
template <typename Kernel>
CGAL::Plane_3<Kernel> bisector(const CGAL::Point_3<Kernel> &p,
const CGAL::Point_3<Kernel> &q);

/*!
constructs the bisector of the two planes `h1` and `h2`.
In the general case, the bisector has a normal vector which has the same
direction as the sum of the normalized normal vectors of the two planes, and
passes through the intersection of `h1` and `h2`.
If `h1` and `h2` are parallel, then the bisector is defined as the
plane which has the same oriented normal vector as `l1`, and which is at
the same distance from `h1` and `h2`.
This function requires that `Kernel::RT` supports the `sqrt()`
operation.
*/
template <typename Kernel>
CGAL::Plane_3<Kernel> bisector(const CGAL::Plane_3<Kernel> &h1,
const CGAL::Plane_3<Kernel> &h2);

/// @}

/// \defgroup centroid_grp CGAL::centroid()
/// \ingroup kernel_global_function
/// \sa \link barycenter `CGAL::barycenter()` \endlink
/// @{

/*!
compute the centroid of the points `p`, `q`, and `r`.
*/
template <typename Kernel>
CGAL::Point_2<Kernel>
centroid( const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r);

/*!
compute the centroid of the points `p`, `q`, `r`, and `s`.
*/
template <typename Kernel>
CGAL::Point_2<Kernel>
centroid( const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r,
const CGAL::Point_2<Kernel>& s);

/*!
compute the centroid of the triangle `t`.
*/
template <typename Kernel>
CGAL::Point_2<Kernel>
centroid( const CGAL::Triangle_2<Kernel>& t);

/*!
compute the centroid of the points `p`, `q`, and `r`.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
centroid( const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r);

/*!
compute the centroid of the points `p`, `q`, `r`, and `s`.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
centroid( const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r,
const CGAL::Point_3<Kernel>& s);

/*!
compute the centroid of the triangle `t`.
*/
template <typename Kernel>template <typename Kernel>
CGAL::Point_3<Kernel>
centroid( const CGAL::Triangle_3<Kernel>& t);

/*!
compute the centroid of the tetrahedron `t`.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
centroid( const CGAL::Tetrahedron_3<Kernel>& t);

/// @}

/// \defgroup circumcenter_grp CGAL::circumcenter()
/// \ingroup kernel_global_function
/// @{

/*!
compute the center of the smallest circle passing through the points `p` and
`q`. Note: this is the same as `CGAL::midpoint(p, q)` but is provided
for homogeneity. 
*/
template <typename Kernel>
CGAL::Point_2<Kernel>
circumcenter( const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q);

/*!
compute the center of the circle passing through the points `p`, `q`, and `r`.
\pre `p`, `q`, and `r` are not collinear.
*/
template <typename Kernel>
CGAL::Point_2<Kernel>
circumcenter( const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r);

/*!
compute the center of the circle passing through the vertices of `t`.
\pre `t` is not degenerate.
*/
template <typename Kernel>
CGAL::Point_2<Kernel>
circumcenter( const CGAL::Triangle_2<Kernel>& t);

/*!
compute the center of the smallest sphere passing through the points `p` and
`q`. Note: this is the same as `CGAL::midpoint(p, q)` but is provided
for homogeneity. 
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
circumcenter( const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q);

/*!
compute the center of the circle passing through the points `p`, `q`, and `r`.
\pre `p`, `q`, and `r` are not collinear.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
circumcenter( const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r);

/*!
compute the center of the circle passing through the vertices of `t`.
\pre `t` is not degenerate.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
circumcenter( const CGAL::Triangle_3<Kernel>& t);

/*!
compute the center of the sphere passing through the points `p`, `q`, `r`, and `s`.
\pre `p`, `q`, `r`, and `s` are not coplanar.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
circumcenter( const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r,
const CGAL::Point_3<Kernel>& s);

/*!
compute the center of the sphere passing through the vertices of `t`.
\pre `t` is not degenerate.
*/
template <typename Kernel>
CGAL::Point_3<Kernel>
circumcenter( const CGAL::Tetrahedron_3<Kernel>& t);

/// @}

/// \defgroup collinear_are_ordered_along_line_grp CGAL::collinear_are_ordered_along_line()
/// \ingroup kernel_global_function
/// \sa `are_ordered_along_line_grp`
/// \sa `are_strictly_ordered_along_line_grp`
/// \sa `collinear_are_strictly_ordered_along_line_grp`
/// @{

/*!
returns `true`, iff `q` lies between `p` 
and `r`. \pre `p, q` and `r` are collinear.
*/
template <typename Kernel>
bool collinear_are_ordered_along_line(const CGAL::Point_2<Kernel> &p,
const CGAL::Point_2<Kernel> &q,
const CGAL::Point_2<Kernel> &r);

/*!
returns `true`, iff `q` lies between `p` 
and `r`. \pre `p, q` and `r` are collinear.
*/
template <typename Kernel>
bool collinear_are_ordered_along_line(const CGAL::Point_3<Kernel> &p,
const CGAL::Point_3<Kernel> &q,
const CGAL::Point_3<Kernel> &r);

/// @}

/// \defgroup collinear_are_strictly_ordered_along_line_grp CGAL::collinear_are_strictly_ordered_along_line()
/// \ingroup kernel_global_function
/// \sa `are_ordered_along_line_grp`
/// \sa `are_strictly_ordered_along_line_grp`
/// \sa `collinear_are_ordered_along_line_grp`
/// @{

/*!
returns `true`, iff `q` lies strictly between 
`p` and `r`. \pre `p, q` and `r` are collinear.
*/
template <typename Kernel>
bool collinear_are_strictly_ordered_along_line(const CGAL::Point_2<Kernel> &p,
const CGAL::Point_2<Kernel> &q,
const CGAL::Point_2<Kernel> &r);

/*!
returns `true`, iff `q` lies strictly between `p` 
and `r`. \pre `p, q` and `r` are collinear.
*/
template <typename Kernel>
bool collinear_are_strictly_ordered_along_line(
const CGAL::Point_3<Kernel> &p,
const CGAL::Point_3<Kernel> &q,
const CGAL::Point_3<Kernel> &r);

/// @}

/// \defgroup collinear_grp CGAL::collinear()
/// \ingroup kernel_global_function
/// \sa `left_turn_grp`
/// \sa `orientation_grp`
/// \sa `right_turn_grp`
/// @{

/*!
returns `true`, iff `p`, `q`, and `r` are collinear.
*/
template <typename Kernel>
bool collinear(const CGAL::Point_2<Kernel> &p, 
const CGAL::Point_2<Kernel> &q, 
const CGAL::Point_2<Kernel> &r);

/*!
returns `true`, iff `p`, `q`, and `r` are collinear.
*/
template <typename Kernel>
bool collinear(const CGAL::Point_3<Kernel> &p,
const CGAL::Point_3<Kernel>&q,
const CGAL::Point_3<Kernel>&r);

/// @}



/// \defgroup compare_dihedral_angle_grp CGAL::compare_dihedral_angle()
/// \ingroup kernel_global_function
/// @{

/*!
compares the dihedral angles \f$ \theta_1\f$ and \f$ \theta_2\f$, where
\f$ \theta_1\f$ is the dihedral angle, in \f$ [0, \pi]\f$, of the tetrahedron
`(a_1, b_1, c_1, d_1)`at the edge `(a_1, b_1)`, and \f$ \theta_2\f$ is
the angle in \f$ [0, \pi]\f$ such that \f$ cos(\theta_2) = cosine\f$.
The result is the same as `compare_dihedral_angle(b1-a1, c1-a1, d1-a1, cosine)`.
\pre `a_1`, `b_1`, `c_1` are not collinear, and `a_1`, `b_1`, `d_1` are not collinear.
*/
template <typename Kernel>
Comparison_result compare_dihedral_angle(const CGAL::Point_3<Kernel>& a1,
const CGAL::Point_3<Kernel>& b1, 
const CGAL::Point_3<Kernel>& c1,
const CGAL::Point_3<Kernel>& d1, 
const Kernel::FT& cosine);

/*!
compares the dihedral angles \f$ \theta_1\f$ and \f$ \theta_2\f$, where
\f$ \theta_i\f$ is the dihedral angle in the tetrahedron `(a_i, b_i,
c_i, d_i)`at the edge `(a_i, b_i)`. These two angles are computed
in \f$ [0, \pi]\f$.
The result is the same as `compare_dihedral_angle(b1-a1, c1-a1, d1-a1, b2-a2, c2-a2, d2-a2)`.
\pre For \f$ i \in\{1,2\}\f$, `a_i`, `b_i`, `c_i` are not collinear, and `a_i`, `b_i`, `d_i` are not collinear.
*/
template <typename Kernel>
Comparison_result compare_dihedral_angle(const CGAL::Point_3<Kernel>& a1,
const CGAL::Point_3<Kernel>& b1, 
const CGAL::Point_3<Kernel>& c1,
const CGAL::Point_3<Kernel>& d1, 
const CGAL::Point_3<Kernel>& a2, 
const CGAL::Point_3<Kernel>& b2, 
const CGAL::Point_3<Kernel>& c2,
const CGAL::Point_3<Kernel>& d2);

/*!
compares the dihedral angles \f$ \theta_1\f$ and \f$ \theta_2\f$, where
\f$ \theta_1\f$ is the dihedral angle, in \f$ [0, \pi]\f$, between the
vectorial planes defined by `(u_1, v_1)` and `(u_1, w_1)`, and
\f$ \theta_2\f$ is the angle in \f$ [0, \pi]\f$ such that \f$ cos(\theta_2) =
cosine\f$.
\pre `u_1` and `v_1` are not collinear, and `u_1` and `w_1` are not collinear.
*/
template <typename Kernel>
Comparison_result
compare_dihedral_angle(const CGAL::Vector_3<Kernel>& u1,
const CGAL::Vector_3<Kernel>& v1, 
const CGAL::Vector_3<Kernel>& w1,
const Kernel::FT& cosine);

/*!
compares the dihedral angles \f$ \theta_1\f$ and \f$ \theta_2\f$, where
\f$ \theta_i\f$ is the dihedral angle between the vectorial planes
defined by  `(u_i, v_i)` and `(u_i, w_i)`. These two angles are
computed in \f$ [0, \pi]\f$.
\pre For \f$ i \in\{1,2\}\f$, `u_i` and `v_i` are not collinear, and `u_i` and `w_i` are not collinear.
*/
template <typename Kernel>
Comparison_result
compare_dihedral_angle(const CGAL::Vector_3<Kernel>& u1,
const CGAL::Vector_3<Kernel>& v1, 
const CGAL::Vector_3<Kernel>& w1,
const CGAL::Vector_3<Kernel>& u2, 
const CGAL::Vector_3<Kernel>& v2, 
const CGAL::Vector_3<Kernel>& w2);

/// @}

/// \defgroup compare_distance_to_point_grp CGAL::compare_distance_to_point()
/// \ingroup kernel_global_function
/// \sa `compare_squared_distance_grp`
/// \sa `compare_signed_distance_to_line_grp`
/// \sa `compare_signed_distance_to_plane_grp`
/// \sa `has_larger_distance_to_point_grp`
/// \sa `has_larger_signed_distance_to_line_grp`
/// \sa `has_larger_signed_distance_to_plane_grp`
/// \sa `has_smaller_distance_to_point_grp`
/// \sa `has_smaller_signed_distance_to_line_grp`
/// \sa `has_smaller_signed_distance_to_plane_grp`
/// @{

/*!
compares the distances of points `q` and
`r` to point `p`.
returns `CGAL::SMALLER`, iff `q` is closer
to `p` than `r`, `CGAL::LARGER`, iff
`r` is closer to `p` than `q`, and
`CGAL::EQUAL` otherwise.
*/
template <typename Kernel>
Comparison_result
compare_distance_to_point(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r);

/*!
compares the distances of points `q` and
`r` to point `p`.
returns `CGAL::SMALLER`, iff `q` is closer
to `p` than `r`, `CGAL::LARGER`, iff
`r` is closer to `p` than `q`, and
`CGAL::EQUAL` otherwise.
*/
template <typename Kernel>
Comparison_result
compare_distance_to_point(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r);

/// @}



/// \defgroup compare_lexicographically_linear_grp CGAL::compare_lexicographically()
/// \ingroup kernel_global_function
/// @{

/*!
Compares the Cartesian coordinates of points `p` and
`q` lexicographically in \f$ xy\f$ order: first 
\f$ x\f$-coordinates are compared, if they are equal, \f$ y\f$-coordinates
are compared. This is the same function as `compare_xy` and exists for compatibility with `Point_d<Kernel>`.
*/
template <typename Kernel>
Comparison_result
compare_lexicographically(const CGAL::Point_2<Kernel>& p, const CGAL::Point_2<Kernel>& q);

/*!
Compares the Cartesian coordinates of points `p` and
`q` lexicographically in \f$ xyz\f$ order: first 
\f$ x\f$-coordinates are compared, if they are equal, \f$ y\f$-coordinates
are compared, and if both \f$ x\f$- and \f$ y\f$- coordinate are equal,
\f$ z\f$-coordinates are compared. This is the same function as `compare_xyz` and exists for compatibility with `Point_d<Kernel>`.
*/
template <typename Kernel>
Comparison_result
compare_lexicographically(const CGAL::Point_3<Kernel>& p, const CGAL::Point_3<Kernel>& q);

/// @}


/// \defgroup compare_signed_distance_to_line_grp CGAL::compare_signed_distance_to_line()
/// \ingroup kernel_global_function
/// \sa `compare_distance_to_point_grp`
/// \sa `compare_signed_distance_to_plane_grp`
/// \sa `has_larger_distance_to_point_grp`
/// \sa `has_larger_signed_distance_to_line_grp`
/// \sa `has_larger_signed_distance_to_plane_grp`
/// \sa `has_smaller_distance_to_point_grp`
/// \sa `has_smaller_signed_distance_to_line_grp`
/// \sa `has_smaller_signed_distance_to_plane_grp`
/// @{

/*!
returns `CGAL::LARGER`
iff the signed distance of `p` and 
`l` is larger than the signed distance of `q`
and `l`, `CGAL::SMALLER`, iff it is smaller,
and `CGAL::EQUAL` iff both are equal.
*/
template <typename Kernel>
Comparison_result
compare_signed_distance_to_line(const CGAL::Line_2<Kernel>& l,
const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q);

/*!
returns `CGAL::LARGER`
iff the signed distance of `r` and 
`l` is larger than the signed distance of `s`
and `l`, `CGAL::SMALLER`, iff it is smaller,
and `CGAL::EQUAL` iff both are equal, where 
`l` is the directed line through `p` and `q`.
*/
template <typename Kernel>
Comparison_result
compare_signed_distance_to_line(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r,
const CGAL::Point_2<Kernel>& s);

/// @}

/// \defgroup compare_signed_distance_to_plane_grp CGAL::compare_signed_distance_to_plane()
/// \ingroup kernel_global_function
/// \sa `compare_distance_to_point_grp`
/// \sa `compare_signed_distance_to_line_grp`
/// \sa `has_larger_distance_to_point_grp`
/// \sa `has_larger_signed_distance_to_line_grp`
/// \sa `has_larger_signed_distance_to_plane_grp`
/// \sa `has_smaller_distance_to_point_grp`
/// \sa `has_smaller_signed_distance_to_line_grp`
/// \sa `has_smaller_signed_distance_to_plane_grp`
/// @{

/*!
returns `CGAL::LARGER`
iff the signed distance of `p` and 
`h` is larger than the signed distance of `q`
and `h`, `CGAL::SMALLER`, iff it is smaller,
and `CGAL::EQUAL` iff both are equal.
*/
template <typename Kernel>
Comparison_result
compare_signed_distance_to_plane(const CGAL::Plane_3<Kernel>& h,
const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q);

/*!
returns `CGAL::LARGER`
iff the signed distance of `s` and 
`h` is larger than the signed distance of `t`
and `h`, `CGAL::SMALLER`, iff it is smaller,
and `CGAL::EQUAL` iff both are equal, where
`h` is the oriented plane through `p`, `q` and
`r`.
*/
template <typename Kernel>
Comparison_result
compare_signed_distance_to_plane(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r,
const CGAL::Point_3<Kernel>& s,
const CGAL::Point_3<Kernel>& t);

/// @}



/// \defgroup compare_slopes_grp CGAL::compare_slopes()
/// \ingroup kernel_global_function
/// @{

/*!
compares the slopes of the lines `l1` and `l2`
*/
template <typename Kernel>
Comparison_result compare_slopes(const CGAL::Line_2<Kernel> &l1,
const CGAL::Line_2<Kernel> &l2);

/*!
compares the slopes of the segments `s1` and `s2`
*/
template <typename Kernel>
Comparison_result compare_slopes(const CGAL::Segment_2<Kernel> &s1,
const CGAL::Segment_2<Kernel> &s2);

/// @}

/// \defgroup compare_squared_distance_grp CGAL::compare_squared_distance()
/// \ingroup kernel_global_function
/// \sa `compare_distance_to_point_grp`
/// \sa `compare_signed_distance_to_line_grp`
/// \sa `compare_signed_distance_to_plane_grp`
/// \sa `has_larger_distance_to_point_grp`
/// \sa `has_larger_signed_distance_to_line_grp`
/// \sa `has_larger_signed_distance_to_plane_grp`
/// \sa `has_smaller_distance_to_point_grp`
/// \sa `has_smaller_signed_distance_to_line_grp`
/// \sa `has_smaller_signed_distance_to_plane_grp`
/// @{

/*!
compares the squared distance of points `p` and
`q` to `d2`.
*/
template <typename Kernel>
Comparison_result
compare_squared_distance(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const typename Kernel::FT& d2);

/*!
compares the squared distance of points `p` and
`q` to `d2`.
*/
template <typename Kernel>
Comparison_result
compare_squared_distance(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const typename Kernel::FT& d2);

/// @}

/// \defgroup compare_squared_radius_grp CGAL::compare_squared_radius()
/// \ingroup kernel_global_function
/// @{

/*!
compares the squared radius of the sphere of radius 0 centered at
`p` to `sr`. 
This returns the opposite sign of `sr`.
*/
template <typename Kernel>
Comparison_result
compare_squared_radius(const CGAL::Point_3<Kernel>& p,
const typename Kernel::FT& sr);

/*!
compares the squared radius of the sphere defined by the
points `p` and `q` to `sr`.
*/
template <typename Kernel>
Comparison_result
compare_squared_radius(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const typename Kernel::FT& sr);

/*!
compares the squared radius of the sphere defined by the
points `p`, `q`, and `r` to `sr`.
*/
template <typename Kernel>
Comparison_result
compare_squared_radius(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r,
const typename Kernel::FT& sr);

/*!
compares the squared radius of the sphere defined by the
points `p`, `q`, `r`, and `r` to `sr`.
*/
template <typename Kernel>
Comparison_result
compare_squared_radius(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r,
const CGAL::Point_3<Kernel>& s,
const typename Kernel::FT& sr);

/// @}

/*!
\defgroup compare_x_grp CGAL::compare_x()
\ingroup kernel_global_function

\details Depending on which \cgal kernel is used,
different versions of this global function are available. This is
described below.

\sa `compare_xy_grp`
\sa `compare_xyz_grp`
\sa `compare_x_at_y_grp`
\sa `compare_y_grp`
\sa `compare_yx_grp`
\sa `compare_y_at_x_grp`
\sa `compare_z_grp`
*/

/*!
\defgroup compare_x_linear_grp CGAL::compare_x() (2D/3D Linear Kernel)
\ingroup compare_x_grp
\anchor figcompare_x
\image html compare1.png 
\image latex compare1.png 
*/
/// @{

/*!
compares the \f$ x\f$-coordinates of `p` and `q`.
*/
template <typename Kernel>
Comparison_result compare_x(const CGAL::Point_2<Kernel> &p, const CGAL::Point_2<Kernel> &q);

/*!
compares the \f$ x\f$-coordinates of `p` and `q`.
*/
template <typename Kernel>
Comparison_result compare_x(const CGAL::Point_3<Kernel> &p, const CGAL::Point_3<Kernel> &q);
/*!
compares the \f$ x\f$-coordinates of `p` and the intersection 
of lines `l1` and `l2`.
See Figure \ref figcompare_x (a).
*/
template <typename Kernel>
Comparison_result compare_x(const CGAL::Point_2<Kernel> &p,
                            const CGAL::Line_2<Kernel> &l1,
                            const CGAL::Line_2<Kernel> &l2);

/*!
compares the \f$ x\f$-coordinates of  the intersection of line `l`
with line `h1` and with line `h2`.

See Figure \ref figcompare_x (b).
*/
template <typename Kernel>
Comparison_result compare_x(const CGAL::Line_2<Kernel> &l,
                            const CGAL::Line_2<Kernel> &h1,
                            const CGAL::Line_2<Kernel> &h2);
/*!
compares the \f$ x\f$-coordinates of the intersection of lines `l1`
and `l2` and  the intersection of lines `h1` and `h2`.

See Figure \ref figcompare_x (c).
*/
template <typename Kernel>
Comparison_result compare_x(const CGAL::Line_2<Kernel> &l1,
                                        const CGAL::Line_2<Kernel> &l2,
                                        const CGAL::Line_2<Kernel> &h1,
                                        const CGAL::Line_2<Kernel> &h2);

/// @}

/*!
\defgroup compare_x_circular_grp CGAL::compare_x() (2D Circular Kernel)
\ingroup compare_x_grp
\details See Chapter \ref Chapter_2D_Circular_Geometry_Kernel "2D Circular Geometry Kernel".

\code
#include <CGAL/global_functions_circular_kernel_2.h>
\endcode
*/
/// @{

/*!
compares the \f$ x\f$-coordinates of `p` and `q`.
*/
template <typename CircularKernel>
Comparison_result 
  compare_x(const CGAL::Circular_arc_point_2<CircularKernel> &p,
            const CGAL::Circular_arc_point_2<CircularKernel> &q);
/*!
compares the \f$ x\f$-coordinates of `p` and `q`.
*/
template <typename CircularKernel>
Comparison_result 
  compare_x(const CGAL::Circular_arc_point_2<CircularKernel> &p,
            const CGAL::Point_2<CircularKernel> &q);

/// @}

/*!
\defgroup compare_x_spherical_grp CGAL::compare_x() (3D Spherical Kernel)
\ingroup compare_x_grp
\details See Chapter \ref Chapter_3D_Spherical_Geometry_Kernel "3D Spherical Geometry Kernel".

\code
#include <CGAL/global_functions_spherical_kernel_3.h>
\endcode
*/
/// @{

/*!
compares the \f$ x\f$-coordinates of `p` and `q`.
*/
template <typename SphericalKernel>
Comparison_result 
compare_x(const CGAL::Circular_arc_point_3<SphericalKernel> &p,
          const CGAL::Circular_arc_point_3<SphericalKernel> &q);

/*!
compares the \f$ x\f$-coordinates of `p` and `q`.
*/
template <typename SphericalKernel>
Comparison_result 
  compare_x(const CGAL::Circular_arc_point_3<SphericalKernel> &p,
            const CGAL::Point_3<SphericalKernel> &q);

/// @}

/*!
\defgroup compare_xy_grp CGAL::compare_xy()
\ingroup kernel_global_function

\details Depending on which \cgal kernel is used, different versions of this
global function are available. 

\sa `compare_xyz_grp`
\sa `compare_x_grp`
\sa `compare_x_at_y_grp`
\sa `compare_y_grp`
\sa `compare_yx_grp`
\sa `compare_y_at_x_grp`
\sa `compare_z_grp`

*/


/*!
\defgroup compare_xy_linear_grp CGAL::compare_xy() (2D/3D Linear Kernel)
\ingroup compare_xy_grp
*/
/// @{

/*!
Compares the Cartesian coordinates of points `p` and
`q` lexicographically in \f$ xy\f$ order: first 
\f$ x\f$-coordinates are compared, if they are equal, \f$ y\f$-coordinates
are compared.
*/
template <typename Kernel>
Comparison_result
compare_xy(const CGAL::Point_2<Kernel>& p, const CGAL::Point_2<Kernel>& q);

/*!
Compares the Cartesian coordinates of points `p` and `q`
lexicographically in \f$ xy\f$ order: first \f$ x\f$-coordinates are
compared, if they are equal, \f$ y\f$-coordinates are compared.

*/
template <typename Kernel>
Comparison_result
compare_xy(const CGAL::Point_3<Kernel>& p, const CGAL::Point_3<Kernel>& q);

/// @}

/*!
\defgroup compare_xy_circular_grp CGAL::compare_xy() (2D Circular Kernel)
\ingroup compare_xy_grp
\details See Chapter \ref Chapter_2D_Circular_Geometry_Kernel "2D Circular Geometry Kernel".

\code
#include <CGAL/global_functions_circular_kernel_2.h>
\endcode
*/
/// @{

/*!
Compares the \f$ x\f$ and \f$ y\f$ Cartesian coordinates of points `p` and
`q` lexicographically.
*/
template <typename CircularKernel>
Comparison_result 
  compare_xy(const CGAL::Circular_arc_point_2<CircularKernel> &p,
            const CGAL::Circular_arc_point_2<CircularKernel> &q);

/*!
Compares the \f$ x\f$ and \f$ y\f$ Cartesian coordinates of points `p` and
`q` lexicographically.
*/
template <typename CircularKernel>
Comparison_result 
compare_xy(const CGAL::Circular_arc_point_2<CircularKernel> &p,
            const CGAL::Point_2<CircularKernel> &q);

/// @}

/*!
\defgroup compare_xy_spherical_grp CGAL::compare_xy() (3D Spherical Kernel)
\ingroup compare_xy_grp
\details See Chapter \ref Chapter_3D_Spherical_Geometry_Kernel "3D Spherical Geometry Kernel".

\code
#include <CGAL/global_functions_spherical_kernel_3.h>
\endcode
*/
/// @{

/*!

Compares the \f$ x\f$ and \f$ y\f$ Cartesian coordinates of points `p` and
`q` lexicographically.
*/
template <typename SphericalKernel>
Comparison_result 
  compare_xy(const CGAL::Circular_arc_point_3<SphericalKernel> &p,
            const CGAL::Circular_arc_point_3<SphericalKernel> &q);
/*!

Compares the \f$ x\f$ and \f$ y\f$ Cartesian coordinates of points `p` and
`q` lexicographically.
*/
template <typename SphericalKernel>
Comparison_result 
  compare_xy(const CGAL::Circular_arc_point_3<SphericalKernel> &p,
            const CGAL::Point_3<SphericalKernel> &q);

/// @}

/*!
\defgroup compare_x_at_y_grp CGAL::compare_x_at_y()
\ingroup kernel_global_function

\anchor figcomparexaty
\image html compare_x_at_y.png
\image latex compare_x_at_y.png

\sa `compare_xy_grp`
\sa `compare_xyz_grp`
\sa `compare_x_grp`
\sa `compare_y_grp`
\sa `compare_yx_grp`
\sa `compare_y_at_x_grp`
\sa `compare_z_grp`
*/
/// @{

/*!
compares the \f$ x\f$-coordinates of `p` and the horizontal projection
of `p` on `h`.
See Figure \ref figcomparexaty (a).
\pre `h` is not horizontal.
*/
template <typename Kernel>
Comparison_result compare_x_at_y(const CGAL::Point_2<Kernel> &p,
const CGAL::Line_2<Kernel> &h);

/*!
This function compares the \f$ x\f$-coordinates of the horizontal projection 
of `p` on `h1` and on `h2`.
See Figure \ref figcomparexaty (b).
\pre `h1` and `h2` are not horizontal.
*/
template <typename Kernel>
Comparison_result compare_x_at_y(const CGAL::Point_2<Kernel> &p,
const CGAL::Line_2<Kernel> &h1,
const CGAL::Line_2<Kernel> &h2);

/*!
Let `p` be the intersection of lines `l1` and `l2`.
This function compares the \f$ x\f$-coordinates of `p` and 
the horizontal projection of `p` on `h`.
See Figure \ref figcomparexaty (c).
\pre `l1` and `l2` intersect and are not horizontal; `h` is not horizontal.
*/
template <typename Kernel>
Comparison_result compare_x_at_y(const CGAL::Line_2<Kernel> &l1,
const CGAL::Line_2<Kernel> &l2,
const CGAL::Line_2<Kernel> &h);

/*!
Let `p` be the intersection of lines `l1` and `l2`. This 
function compares the \f$ x\f$-coordinates of the horizontal projection of 
`p` on `h1` and on `h2`
See Figure \ref figcomparexaty (d).
\pre `l1` and `l2` intersect and are not horizontal; `h1` and `h2` are not horizontal.
*/
template <typename Kernel>
Comparison_result compare_x_at_y(const CGAL::Line_2<Kernel> &l1,
const CGAL::Line_2<Kernel> &l2,
const CGAL::Line_2<Kernel> &h1,
const CGAL::Line_2<Kernel> &h2);

/// @}

/*!
  \defgroup compare_y_at_x_grp CGAL::compare_y_at_x()
  \ingroup kernel_global_function

  \anchor figcompareyatx
  \image html compare2.png
  \image latex compare2.png

  \sa `compare_xy_grp`
  \sa `compare_xyz_grp`
  \sa `compare_x_grp`
  \sa `compare_y_grp`
  \sa `compare_yx_grp`
  \sa `compare_x_at_y_grp`
  \sa `compare_z_grp`
*/
/// @{

/*!
  compares the \f$ y\f$-coordinates of  `p` and the vertical projection
  of `p` on `h`.
  See Figure \ref figcompareyatx (d).

  \pre `h` is not vertical.
*/
template <typename Kernel>
Comparison_result compare_y_at_x(const CGAL::Point_2<Kernel> &p,
                                 const CGAL::Line_2<Kernel> &h);

/*!
  compares the \f$ y\f$-coordinates of the vertical projection 
  of `p` on `h1` and on `h2`.
  See Figure \ref figcompareyatx (e).

  \pre `h1` and `h2` are not vertical.
*/
template <typename Kernel>
Comparison_result compare_y_at_x(const CGAL::Point_2<Kernel> &p,
                                 const CGAL::Line_2<Kernel> &h1,
                                 const CGAL::Line_2<Kernel> &h2);


/*!
  Let `p` be the `intersection` of lines `l1` and `l2`.
  This function compares the \f$ y\f$-coordinates of `p` and 
  the vertical projection of `p` on `h`
  See Figure \ref figcompareyatx (f).

  \pre `l1`, `l2` intersect and `h` is not vertical.
*/
template <typename Kernel>
Comparison_result compare_y_at_x(const CGAL::Line_2<Kernel> &l1,
                                 const CGAL::Line_2<Kernel> &l2,
                                 const CGAL::Line_2<Kernel> &h);

/*!
  Let `p` be the `intersection` of lines `l1` and `l2`. This function 
  compares the \f$ y\f$-coordinates of the vertical projection of `p` on 
  `h1` and on `h2`.
  See Figure \ref figcompareyatx (g).
  \pre `l1` and `l2` intersect; `h1` and  `h2` are not vertical.
*/
template <typename Kernel>
Comparison_result compare_y_at_x(const CGAL::Line_2<Kernel> &l1,
                                 const CGAL::Line_2<Kernel> &l2,
                                 const CGAL::Line_2<Kernel> &h1,
                                 const CGAL::Line_2<Kernel> &h2);

/*!
  compares the \f$ y\f$-coordinates of `p` and the vertical projection
  of `p` on `s`.  If `s` is vertical, then return
  `CGAL::EQUAL` when `p` lies on `s`, `CGAL::SMALLER` when `p` lies
  under {s}, and `CGAL::LARGER` otherwise.
  \pre `p` is within the x range of `s`.
*/
template <typename Kernel>
Comparison_result compare_y_at_x(const CGAL::Point_2<Kernel> &p,
                                 const CGAL::Segment_2<Kernel> &s);

/*!
  compares the \f$ y\f$-coordinates of the vertical projection 
  of `p` on `s1` and on `s2`.  If `s1` or `s2`
  is vertical, then return `CGAL::EQUAL` if they intersect, otherwise return
  `CGAL::SMALLER` if `s1` lies below `s2`, and return `CGAL::LARGER`
  otherwise.
  \pre `p` is within the x range of `s1` and `s2`.
*/
template <typename Kernel>
Comparison_result compare_y_at_x(const CGAL::Point_2<Kernel> &p,
                                 const CGAL::Segment_2<Kernel> &s1,
                                 const CGAL::Segment_2<Kernel> &s2);

/*!
  \name With the 2D Circular Kernel
  See \ref Chapter_2D_Circular_Geometry_Kernel "2D Circular Geometry Kernel".

  \code 
  #include <CGAL/global_functions_circular_kernel_2.h>
  \endcode
*/
/// @{

/// Same as above, for a point and a circular arc.
template <typename CircularKernel>
Comparison_result 
compare_y_at_x(const CGAL::Circular_arc_point_2<CircularKernel> &p, 
               const CGAL::Circular_arc_2<CircularKernel> &a);

/// Same as above, for a point and a line segment.
template <typename CircularKernel>
Comparison_result 
compare_y_at_x(const CGAL::Circular_arc_point_2<CircularKernel> &p, 
               const CGAL::Line_arc_2<CircularKernel> &a);


/// @}

/// @}


/*!
\defgroup compare_y_grp CGAL::compare_y()
\ingroup kernel_global_function

\details Depending on which \cgal kernel is used, different versions of this
global function are available.

\sa `compare_xy_grp`
\sa `compare_xyz_grp`
\sa `compare_x_grp`
\sa `compare_x_at_y_grp`
\sa `compare_yx_grp`
\sa `compare_y_at_x_grp`
\sa `compare_z_grp`
*/

/*!
\defgroup compary_y_linear_grp CGAL::compare_y() (2D/3D Linear Kernel)
\ingroup compare_y_grp
\details See Chapter \ref chapterkernel23 "2D and 3D Geometry Kernel"

\anchor figcompare13
\image html compare1.png
\image latex compare1.png

*/
/// @{
/*!
  compares Cartesian \f$ y\f$-coordinates of `p` and `q`.
*/
template <typename Kernel>
Comparison_result compare_y(const CGAL::Point_2<Kernel> &p,
                            const CGAL::Point_2<Kernel> &q);
/*!
  compares Cartesian \f$ y\f$-coordinates of `p` and `q`.
*/
template <typename Kernel>
Comparison_result compare_y(const CGAL::Point_3<Kernel> &p,
                            const CGAL::Point_3<Kernel> &q);

/*!
  compares the \f$ y\f$-coordinates of `p` and the intersection of lines
  `l1` and `l2`.
  See Figure \ref figcompare13 (a).
*/
template <typename Kernel>
Comparison_result compare_y(const CGAL::Point_2<Kernel> &p,
                            const CGAL::Line_2<Kernel> &l1,
                            const CGAL::Line_2<Kernel> &l2);
/*!

  compares the \f$ y\f$-coordinates of the intersection of line `l`
  with line `h1` and with line `h2`.
  See Figure \ref figcompare13 (b).
*/
template <typename Kernel>
Comparison_result compare_y(const CGAL::Line_2<Kernel> &l,
                            const CGAL::Line_2<Kernel> &h1,
                            const CGAL::Line_2<Kernel> &h2);
/*!
  compares the \f$ y\f$-coordinates of the intersection of lines `l1`
  and `l2` and  the intersection of lines `h1` and `h2`.
  See Figure \ref figcompare13 (c).
*/
template <typename Kernel>
Comparison_result compare_y(const CGAL::Line_2<Kernel> &l1,
                            const CGAL::Line_2<Kernel> &l2,
                            const CGAL::Line_2<Kernel> &h1,
                            const CGAL::Line_2<Kernel> &h2);

/// @}

/*!
\defgroup compare_y_circular_grp CGAL::compare_y() (2D Circular Kernel)
\ingroup compare_y_grp
\details See Chapter \ref Chapter_2D_Circular_Geometry_Kernel "2D Circular Geometry Kernel".

\code
#include <CGAL/global_functions_circular_kernel_2.h>
\endcode
*/
/// @{
/*!
  compares the \f$ y\f$-coordinates of `p` and `q`.
*/
template <typename CircularKernel>
Comparison_result 
compare_y(const CGAL::Circular_arc_point_2<CircularKernel> &p,
          const CGAL::Circular_arc_point_2<CircularKernel> &q);

/*!
  compares the \f$ y\f$-coordinates of `p` and `q`.
*/template <typename CircularKernel>
Comparison_result 
compare_y(const CGAL::Circular_arc_point_2<CircularKernel> &p,
          const CGAL::Point_2<CircularKernel> &q);

/// @}

/*!
\defgroup compare_y_spherical_grp CGAL::compare_y() (3D Spherical Kernel)
\ingroup compare_y_grp
\details See Chapter \ref Chapter_3D_Spherical_Geometry_Kernel "3D Spherical Geometry Kernel".

\code
#include <CGAL/global_functions_circular_kernel_3.h>
\endcode
*/
/// @{
/*!
compares the \f$ y\f$-coordinates of `p` and `q`.
*/
template <typename SphericalKernel>
Comparison_result 
  compare_y(const CGAL::Circular_arc_point_3<SphericalKernel> &p,
            const CGAL::Circular_arc_point_3<SphericalKernel> &q);
/*!
compares the \f$ y\f$-coordinates of `p` and `q`.
*/
template <typename SphericalKernel>
Comparison_result 
  compare_y(const CGAL::Circular_arc_point_3<SphericalKernel> &p,
            const CGAL::Point_3<SphericalKernel> &q);
/// @}


/*!
\defgroup compare_xyz_grp CGAL::compare_xyz()
\ingroup kernel_global_function

\details Depending on which \cgal kernel is used, different versions of this
global function are available.

\sa `compare_xy_grp`
\sa `compare_x_grp`
\sa `compare_x_at_y_grp`
\sa `compare_y_grp`
\sa `compare_yx_grp`
\sa `compare_y_at_x_grp`
\sa `compare_z_grp`

*/

/*!
\defgroup compare_xyz_linear_grp CGAL::compare_xyz() (2D/3D Linear Kernel)
\ingroup compare_xyz_grp
*/
/// @{

/*!
Compares the Cartesian coordinates of points `p` and
`q` lexicographically in \f$ xyz\f$ order: first 
\f$ x\f$-coordinates are compared, if they are equal, \f$ y\f$-coordinates
are compared, and if both \f$ x\f$- and \f$ y\f$- coordinate are equal,
\f$ z\f$-coordinates are compared.
*/
template <typename Kernel>
Comparison_result
compare_xyz(const CGAL::Point_3<Kernel>& p, const CGAL::Point_3<Kernel>& q);

/// @}

/*!
\defgroup compare_xyz_spherical_grp CGAL::compare_xyz() (3D Spherical Kernel)
\ingroup compare_xyz_grp
\details See Chapter \ref Chapter_3D_Spherical_Geometry_Kernel "3D Spherical Geometry Kernel"

\code
#include <CGAL/global_functions_circular_kernel_3.h>
\endcode
*/
/// @{

/*! Compares the Cartesian coordinates of points `p` and `q` lexicographically.
*/
template <typename SphericalKernel>
Comparison_result 
compare_xyz(const CGAL::Circular_arc_point_3<SphericalKernel> &p,
const CGAL::Circular_arc_point_3<SphericalKernel> &q);

/*!
Compares the Cartesian coordinates of points `p` and `q` lexicographically.
*/
template <typename SphericalKernel>
Comparison_result 
compare_xyz(const CGAL::Circular_arc_point_3<SphericalKernel> &p,
const CGAL::Point_3<SphericalKernel> &q);

/// @}


/*!
\defgroup compare_z_grp CGAL::compare_z()
\ingroup kernel_global_function

\details Depending on which \cgal kernel is used, 
different versions of this global function are available. This is 
described below. 

\sa `compare_xy_grp`
\sa `compare_xyz_grp`
\sa `compare_x_grp`
\sa `compare_x_at_y_grp`
\sa `compare_y_grp`
\sa `compare_yx_grp`
\sa `compare_y_at_x_grp`
*/

/*!
\defgroup compare_z_linear_grp CGAL::compare_z() (2D/3D Linear Kernel)
\ingroup compare_z_grp
*/
/// @{

/*!
compares the \f$ z\f$-coordinates of `p` and `q`.
*/
template <typename Kernel>
Comparison_result compare_z(const CGAL::Point_3<Kernel> &p, const CGAL::Point_3<Kernel> &q);

/// @}

/*!
\defgroup compare_z_spherical_grp CGAL::compare_z() (3D Spherical Kernel)
\ingroup compare_z_grp

\details See Chapter \ref Chapter_3D_Spherical_Geometry_Kernel "3D Spherical Geometry Kernel"

\code
#include <CGAL/global_functions_circular_kernel_3.h>
\endcode

*/
/*!



*/
/// @{

/*!
compares the \f$ z\f$-coordinates of `p` and `q`.
*/
template <typename SphericalKernel>
Comparison_result 
compare_z(const CGAL::Circular_arc_point_3<SphericalKernel> &p, const CGAL::Circular_arc_point_3<SphericalKernel> &q);

/*!
compares the \f$ z\f$-coordinates of `p` and `q`.
*/
template <typename SphericalKernel>
Comparison_result 
compare_z(const CGAL::Circular_arc_point_3<SphericalKernel> &p, const CGAL::Point_3<SphericalKernel> &q);

/// @}

/// \defgroup compare_yx_grp CGAL::compare_yx()
/// \ingroup kernel_global_function
/// \sa `compare_xy_grp`
/// \sa `compare_xyz_grp`
/// \sa `compare_x_grp`
/// \sa `compare_x_at_y_grp`
/// \sa `compare_y_grp`
/// \sa `compare_y_at_x_grp`
/// \sa `compare_z_grp`
/// @{

/*!
Compares the Cartesian coordinates of points `p` and
`q` lexicographically in \f$ yx\f$ order: first 
\f$ y\f$-coordinates are compared, if they are equal, \f$ x\f$-coordinates
are compared.
*/
template <typename Kernel>
Comparison_result
compare_yx(const CGAL::Point_2<Kernel>& p, const CGAL::Point_2<Kernel>& q);

/// @}


/// \defgroup coplanar_grp CGAL::coplanar()
/// \ingroup kernel_global_function
/// \sa `coplanar_orientation_grp`
/// \sa `coplanar_side_of_bounded_circle_grp`
/// @{

/*!
returns `true`, if `p`, `q`, `r`, and `s` are coplanar.
*/
template <typename Kernel>
bool coplanar(const CGAL::Point_3<Kernel> &p,
const CGAL::Point_3<Kernel>&q,
const CGAL::Point_3<Kernel>&r,
const CGAL::Point_3<Kernel>&s);

/// @}

/// \defgroup coplanar_orientation_grp CGAL::coplanar_orientation()
/// \ingroup kernel_global_function
/// \sa `coplanar_grp`
/// \sa `coplanar_side_of_bounded_circle_grp`
/// \sa `orientation_grp`
/// @{

/*!
Let `p` be the plane defined by the points `p`, `q`,
and `r`. Note that the order defines the orientation of
`p`. The function computes the orientation of points `p`, 
`q`, and `s` in `p`: Iff `p`, `q`, `s` are
collinear, `CGAL::COLLINEAR` is returned. Iff `p` and the plane 
defined by `p`, `q`, and `s` have the same orientation, 
`CGAL::POSITIVE` is returned; otherwise `CGAL::NEGATIVE` is returned. 
\pre `p`, `q`, `r`, and `s` are coplanar and `p`, `q`, and `r` are not collinear.
*/
template <typename Kernel>
Orientation coplanar_orientation(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r,
const CGAL::Point_3<Kernel>& s);

/*!
If `p,q,r` are collinear, then `CGAL::COLLINEAR` is returned.
If not, then `p,q,r` define a plane `p`. The return value in this case is
either `CGAL::POSITIVE` or `CGAL::NEGATIVE`, but we don't specify it explicitly.
However, we guarantee that all calls to this predicate over 3 points in `p`
will return a coherent orientation if considered a 2D orientation in `p`.
*/
template <typename Kernel>
Orientation coplanar_orientation(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r);

/// @}



/// \defgroup coplanar_side_of_bounded_circle_grp CGAL::coplanar_side_of_bounded_circle()
/// \ingroup kernel_global_function
/// \sa `coplanar_orientation_grp`
/// \sa `side_of_bounded_circle_grp`
/// @{

/*!
returns the bounded side of the circle defined
by `p`, `q`, and `r` on which `s` lies.
\pre `p`, `q`, `r`, and `s` are coplanar and `p`, `q`, and `r` are not collinear.
*/
template <typename Kernel>
Bounded_side coplanar_side_of_bounded_circle(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r,
const CGAL::Point_3<Kernel>& s);

/// @}


/// \defgroup cross_product_grp CGAL::cross_product()
/// \ingroup kernel_global_function
/// @{

/*!
returns the cross product of `u` and `v`.
*/
template <typename Kernel>
CGAL::Vector_3<Kernel> cross_product( const CGAL::Vector_3<Kernel>& u, 
const CGAL::Vector_3<Kernel>& v);

/// @}

/// \defgroup determinant_grp CGAL::determinant()
/// \ingroup kernel_global_function
/// \sa `orientation_grp`
/// \sa `collinear_grp`
/// \sa `left_turn_grp`
/// \sa `right_turn_grp`
/// @{

/*!
returns the determinant of `v` and `w`.
*/
template <typename Kernel>
Kernel::FT determinant(const CGAL::Vector_2<Kernel>& v,
const CGAL::Vector_2<Kernel>& w);

/*!
returns the determinant of `u`, `v` and `w`.
*/
template <typename Kernel>
Kernel::FT determinant(const CGAL::Vector_3<Kernel>& u,
const CGAL::Vector_3<Kernel>& v,
const CGAL::Vector_3<Kernel>& w);

/// @}

// This is there to keep the global functions in alphabetical order
// instead of processing order.

/// \defgroup do_intersect_grp CGAL::do_intersect()
/// \ingroup kernel_global_function
/// \defgroup do_intersect_linear_grp CGAL::do_intersect() (2D/3D Linear Kernel)
/// \ingroup do_intersect_grp
/// \defgroup do_intersect_circular_grp CGAL::do_intersect() (2D Circular Kernel)
/// \ingroup do_intersect_grp
/// \defgroup do_intersect_spherical_grp CGAL::do_intersect() (3D Spherical Kernel)
/// \ingroup do_intersect_grp

/// \defgroup equidistant_line_grp CGAL::equidistant_line()
/// \ingroup kernel_global_function
/// @{

/*!
constructs the line which is at the same distance from the three points
`p`, `q` and `r`.
\pre `p`, `q` and `r` are not collinear.
*/
template <typename Kernel>
CGAL::Line_3<Kernel> equidistant_line(const CGAL::Point_3<Kernel> &p,
const CGAL::Point_3<Kernel> &q,
const CGAL::Point_3<Kernel> &r);

/// @}

/// \defgroup has_larger_distance_to_point_grp CGAL::has_larger_distance_to_point()
/// \ingroup kernel_global_function
/// \sa `compare_distance_to_point_grp`
/// \sa `compare_signed_distance_to_line_grp`
/// \sa `compare_signed_distance_to_plane_grp`
/// \sa `has_larger_signed_distance_to_line_grp`
/// \sa `has_larger_signed_distance_to_plane_grp`
/// \sa `has_smaller_distance_to_point_grp`
/// \sa `has_smaller_signed_distance_to_line_grp`
/// \sa `has_smaller_signed_distance_to_plane_grp`
/// @{

/*!
returns `true` iff the distance between `q`
and `p` is larger than the distance between `r`
and `p`.
*/
template <typename Kernel>
bool
has_larger_distance_to_point(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r);

/*!
returns `true` iff the distance between `q`
and `p` is larger than the distance between `r`
and `p`.
*/
template <typename Kernel>
bool
has_larger_distance_to_point(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r);

/// @}

/// \defgroup has_larger_signed_distance_to_line_grp CGAL::has_larger_signed_distance_to_line()
/// \ingroup kernel_global_function
/// \sa `compare_distance_to_point_grp`
/// \sa `compare_signed_distance_to_line_grp`
/// \sa `compare_signed_distance_to_plane_grp`
/// \sa `has_larger_distance_to_point_grp`
/// \sa `has_larger_signed_distance_to_line_grp`
/// \sa `has_larger_signed_distance_to_plane_grp`
/// \sa `has_smaller_distance_to_point_grp`
/// \sa `has_smaller_signed_distance_to_line_grp`
/// \sa `has_smaller_signed_distance_to_plane_grp`
/// @{

/*!
returns `true` iff the signed distance of `p`
and `l` is larger than the signed distance of 
`q` and `l`.
*/
template <typename Kernel>
bool
has_larger_signed_distance_to_line(const CGAL::Line_2<Kernel>& l,
const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q);

/*!
returns `true` iff the signed distance of `r`
and `l` is larger than the signed distance of 
`s` and `l`, where `l` is the directed line
through points `p` and `q`.
*/
template <typename Kernel>
bool
has_larger_signed_distance_to_line(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r,
const CGAL::Point_2<Kernel>& s);

/// @}


/// \defgroup has_larger_signed_distance_to_plane_grp CGAL::has_larger_signed_distance_to_plane()
/// \ingroup kernel_global_function
/// \sa `compare_distance_to_point_grp`
/// \sa `compare_signed_distance_to_line_grp`
/// \sa `compare_signed_distance_to_plane_grp`
/// \sa `has_larger_distance_to_point_grp`
/// \sa `has_larger_signed_distance_to_line_grp`
/// \sa `has_smaller_distance_to_point_grp`
/// \sa `has_smaller_signed_distance_to_line_grp`
/// \sa `has_smaller_signed_distance_to_plane_grp`
/// @{

/*!
returns `true` iff the signed distance of `p`
and `h` is larger than the signed distance of 
`q` and `h`.
*/
template <typename Kernel>
bool
has_larger_signed_distance_to_plane(const CGAL::Plane_3<Kernel>& h,
const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q);

/*!
returns `true` iff the signed distance of `s`
and `h` is larger than the signed distance of 
`t` and `h`, where `h` is the oriented
plane through `p`, `q` and `r`.
*/
template <typename Kernel>
bool
has_larger_signed_distance_to_plane(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r,
const CGAL::Point_3<Kernel>& s,
const CGAL::Point_3<Kernel>& t);

/// @}

/// \defgroup has_smaller_distance_to_point_grp CGAL::has_smaller_distance_to_point()
/// \ingroup kernel_global_function
/// \sa `compare_distance_to_point_grp`
/// \sa `compare_signed_distance_to_line_grp`
/// \sa `compare_signed_distance_to_plane_grp`
/// \sa `has_larger_distance_to_point_grp`
/// \sa `has_larger_signed_distance_to_line_grp`
/// \sa `has_larger_signed_distance_to_plane_grp`
/// \sa `has_smaller_signed_distance_to_line_grp`
/// \sa `has_smaller_signed_distance_to_plane_grp`
/// @{

/*!
returns `true` iff the distance between `q`
and `p` is smaller than the distance between `r`
and `p`.
*/
template <typename Kernel>
bool
has_smaller_distance_to_point(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r);

/*!
returns `true` iff the distance between `q`
and `p` is smaller than the distance between `r`
and `p`.
*/
template <typename Kernel>
bool
has_smaller_distance_to_point(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r);

/// @}

/// \defgroup has_smaller_signed_distance_to_line_grp CGAL::has_smaller_signed_distance_to_line()
/// \ingroup kernel_global_function
/// \sa `compare_distance_to_point_grp`
/// \sa `compare_signed_distance_to_line_grp`
/// \sa `compare_signed_distance_to_plane_grp`
/// \sa `has_larger_distance_to_point_grp`
/// \sa `has_larger_signed_distance_to_line_grp`
/// \sa `has_larger_signed_distance_to_plane_grp`
/// \sa `has_smaller_distance_to_point_grp`
/// \sa `has_smaller_signed_distance_to_plane_grp`
/// @{

/*!
returns `true` iff the signed distance of `p`
and `l` is smaller than the signed distance of 
`q` and `l`.
*/
template <typename Kernel>
bool
has_smaller_signed_distance_to_line(const CGAL::Line_2<Kernel>& l,
const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q);

/*!
returns `true` iff the signed distance of `r`
and `l` is smaller than the signed distance of 
`s` and `l`, where `l` is the 
oriented line through `p` and `q`.
*/
template <typename Kernel>
bool
has_smaller_signed_distance_to_line(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r,
const CGAL::Point_2<Kernel>& s);

/// @}

/// \defgroup has_smaller_signed_distance_to_plane_grp CGAL::has_smaller_signed_distance_to_plane()
/// \ingroup kernel_global_function
/// \sa `compare_distance_to_point_grp`
/// \sa `compare_signed_distance_to_line_grp`
/// \sa `compare_signed_distance_to_plane_grp`
/// \sa `has_larger_distance_to_point_grp`
/// \sa `has_larger_signed_distance_to_line_grp`
/// \sa `has_larger_signed_distance_to_plane_grp`
/// \sa `has_smaller_distance_to_point_grp`
/// \sa `has_smaller_signed_distance_to_line_grp`
/// @{

/*!
returns `true` iff the signed distance of `p`
and `h` is smaller than the signed distance of 
`q` and `h`.
*/
template <typename Kernel>
bool
has_smaller_signed_distance_to_plane(const CGAL::Plane_3<Kernel>& h,
const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q);

/*!
returns `true` iff the signed distance of `p`
and `h` is smaller than the signed distance of 
`q` and `h`, where `h` is the oriented
plane through `p`, `q` and `r`.
*/
template <typename Kernel>
bool
has_smaller_signed_distance_to_plane(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r,
const CGAL::Point_3<Kernel>& s,
const CGAL::Point_3<Kernel>& t);

/// @}


// Same reason as in defgroup do_intersect.

/// \defgroup intersection_grp CGAL::intersection()
/// \ingroup kernel_global_function
/// \defgroup intersection_linear_grp CGAL::intersection() (2D/3D Linear Kernel)
/// \ingroup intersection_grp
/// \defgroup intersection_circular_grp CGAL::intersection() (2D Circular Kernel)
/// \ingroup intersection_grp
/// \defgroup intersection_spherical_grp CGAL::intersection() (3D Spherical Kernel)
/// \ingroup intersection_grp

/// \defgroup left_turn_grp CGAL::left_turn()
/// \ingroup kernel_global_function
/// \sa `collinear_grp`
/// \sa `orientation_grp`
/// \sa `right_turn_grp`

/// @{

/*!
returns `true` iff `p`, `q`, and `r` form a left turn.
*/
template <typename Kernel>
bool left_turn(const CGAL::Point_2<Kernel> &p,
const CGAL::Point_2<Kernel> &q,
const CGAL::Point_2<Kernel> &r);

/// @}



/// \defgroup lexicographically_xy_larger_grp CGAL::lexicographically_xy_larger()
/// \ingroup kernel_global_function
/// \sa `compare_xy_grp`
/// \sa `lexicographically_xy_larger_or_equal_grp`
/// \sa `lexicographically_xy_smaller_grp`
/// \sa `lexicographically_xy_smaller_or_equal_grp`
/// @{

/*!
returns `true` iff `p` is lexicographically larger
than `q` with respect to \f$ xy\f$ order.
*/
template <typename Kernel>
bool
lexicographically_xy_larger(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q);

/// @}



/// \defgroup lexicographically_xy_larger_or_equal_grp CGAL::lexicographically_xy_larger_or_equal()
/// \ingroup kernel_global_function
/// \sa `compare_xy_grp`
/// \sa `lexicographically_xy_larger_grp`
/// \sa `lexicographically_xy_smaller_grp`
/// \sa `lexicographically_xy_smaller_or_equal_grp`

/// @{

/*!
returns `true` iff `p` is lexicographically not smaller
than `q` with respect to \f$ xy\f$ order.
*/
template <typename Kernel>
bool
lexicographically_xy_larger_or_equal(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q);

/// @}

/// \defgroup lexicographically_xy_smaller_grp CGAL::lexicographically_xy_smaller()
/// \ingroup kernel_global_function
/// \sa `compare_xy_grp`
/// \sa `lexicographically_xy_larger_grp`
/// \sa `lexicographically_xy_larger_or_equal_grp`
/// \sa `lexicographically_xy_smaller_or_equal_grp`

/// @{

/*!
returns `true` iff `p` is lexicographically smaller
than `q` with respect to \f$ xy\f$ order.
*/
template <typename Kernel>
bool
lexicographically_xy_smaller(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q);

/// @}


/// \defgroup lexicographically_xy_smaller_or_equal_grp CGAL::lexicographically_xy_smaller_or_equal()
/// \ingroup kernel_global_function
/// \sa `compare_xy_grp`
/// \sa `lexicographically_xy_larger_grp`
/// \sa `lexicographically_xy_larger_or_equal_grp`
/// \sa `lexicographically_xy_smaller_grp`

/// @{

/*!
returns `true` iff `p` is lexicographically not larger
than `q` with respect to \f$ xy\f$ order.
*/
template <typename Kernel>
bool 
lexicographically_xy_smaller_or_equal(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q);

/// @}

/// \defgroup lexicographically_xyz_smaller_grp CGAL::lexicographically_xyz_smaller()
/// \ingroup kernel_global_function
/// \sa `compare_xyz_grp`
/// \sa `lexicographically_xyz_smaller_or_equal_grp`

/// @{

/*!
returns `true` iff `p` is lexicographically smaller
than `q` with respect to \f$ xyz\f$ order.
*/
template <typename Kernel>
bool
lexicographically_xyz_smaller(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q);

/// @}

/// \defgroup lexicographically_xyz_smaller_or_equal_grp CGAL::lexicographically_xyz_smaller_or_equal()
/// \ingroup kernel_global_function
/// \sa `compare_xyz_grp`
/// \sa `lexicographically_xyz_smaller_grp`

/// @{

/*!
returns `true` iff `p` is lexicographically not larger
than `q` with respect to \f$ xyz\f$ order.
*/
template <typename Kernel>
bool 
lexicographically_xyz_smaller_or_equal(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q);

/// @}

/// \defgroup max_vertex_grp CGAL::max_vertex()
/// \ingroup kernel_global_function
/// @{

/*!
computes the vertex with the lexicographically largest coordinates of the iso rectangle `ir`.
*/
template <typename Kernel>
CGAL::Point_2<Kernel> max_vertex( const CGAL::Iso_rectangle_2<Kernel>& ir );

/*!
computes the vertex with the lexicographically largest coordinates of the iso cuboid `ic`.
*/
template <typename Kernel>
CGAL::Point_3<Kernel> max_vertex( const CGAL::Iso_cuboid_3<Kernel>& ic );

/// @}

/// \defgroup midpoint_grp CGAL::midpoint()
/// \ingroup kernel_global_function
/// @{

/*!
computes the midpoint of the segment `pq`.
*/
template <typename Kernel>
CGAL::Point_2<Kernel> midpoint( const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q );

/*!
computes the midpoint of the segment `pq`.
*/
template <typename Kernel>
CGAL::Point_3<Kernel> midpoint( const CGAL::Point_3<Kernel>& p, const CGAL::Point_3<Kernel>& q );

/// @}

/// \defgroup min_vertex_grp CGAL::min_vertex()
/// \ingroup kernel_global_function
/// @{

/*!
computes the vertex with the lexicographically smallest coordinates of the iso rectangle `ir`.
*/
template <typename Kernel>
CGAL::Point_2<Kernel> min_vertex( const CGAL::Iso_rectangle_2<Kernel>& ir );

/*!
computes the vertex with the lexicographically smallest coordinates of the iso cuboid `ic`.
*/
template <typename Kernel>
CGAL::Point_3<Kernel> min_vertex( const CGAL::Iso_cuboid_3<Kernel>& ic );

/// @}

/// \defgroup normal_grp CGAL::normal()
/// \ingroup kernel_global_function
/// @{

/*!
computes the normal vector for the vectors `q-p` and `r-p`.
\pre The points `p`, `q`, and `r` must not be collinear.
*/
template <typename Kernel>
CGAL::Vector_3<Kernel> normal( const CGAL::Point_3<Kernel>& p, const CGAL::Point_3<Kernel>& q, const CGAL::Point_3<Kernel>& r );

/// @}

/// \defgroup orientation_grp CGAL::orientation()
/// \ingroup kernel_global_function
/// \sa `collinear_grp`
/// \sa `left_turn_grp`
/// \sa `right_turn_grp`

/// @{

/*!
returns `CGAL::LEFT_TURN`, if `r` lies to the left of the oriented 
line `l` defined by `p` and `q`, returns `CGAL::RIGHT_TURN` if `r` 
lies to the right of `l`, and returns `CGAL::COLLINEAR` if `r` lies
on `l`.
*/
template <typename Kernel>
Orientation orientation(const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r);

/*!
returns `CGAL::LEFT_TURN` if `u` and `v` form a left turn,
returns `CGAL::RIGHT_TURN` if `u` and `v` form a right turn,
and returns `CGAL::COLLINEAR` if `u` and `v` are collinear.
*/
template <typename Kernel>
Orientation orientation(const CGAL::Vector_2<Kernel>& u,
const CGAL::Vector_2<Kernel>& v);

/*!
returns `CGAL::POSITIVE`, if `s` lies on the positive side of the oriented 
plane `h` defined by `p`, `q`, and `r`, returns `CGAL::NEGATIVE` if `s` 
lies on the negative side of `h`, and returns `CGAL::COPLANAR` if `s` lies
on `h`.
*/
template <typename Kernel>
Orientation orientation(const CGAL::Point_3<Kernel> &p,
const CGAL::Point_3<Kernel>&q,
const CGAL::Point_3<Kernel>&r,
const CGAL::Point_3<Kernel>&s);

/*!
returns `CGAL::NEGATIVE` if `u`, `v` and `w` are negatively oriented,
`CGAL::POSITIVE` if `u`, `v` and `w` are positively oriented,
and `CGAL::COPLANAR` if `u`, `v` and `w` are coplanar.
*/
template <typename Kernel>
Orientation orientation(const CGAL::Vector_3<Kernel> &u,
const CGAL::Vector_3<Kernel> &v,
const CGAL::Vector_3<Kernel> &w);

/// @}



/// \defgroup orthogonal_vector_grp CGAL::orthogonal_vector()
/// \ingroup kernel_global_function

/// @{

/*!
computes an orthogonal vector of the plane `p`, which is directed to 
the positive side of this plane.
*/
template <typename Kernel>
CGAL::Vector_3<Kernel> orthogonal_vector( const CGAL::Plane_3<Kernel>& p);

/*!
computes an orthogonal vector of the plane defined by `p`, `q` and `r`,
which is directed to the positive side of this plane.
*/
template <typename Kernel>
CGAL::Vector_3<Kernel> orthogonal_vector( const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r );

/// @}


/// \defgroup parallel_grp CGAL::parallel()
/// \ingroup kernel_global_function

/// @{

/*!
returns `true`, if `l1` and `l2` are parallel or if one
of those (or both) is degenerate.
*/
template <typename Kernel>
bool parallel(const CGAL::Line_2<Kernel>& l1,
const CGAL::Line_2<Kernel>& l2);

/*!
returns `true`, if `r1` and `r2` are parallel or if one
of those (or both) is degenerate.
*/
template <typename Kernel>
bool parallel(const CGAL::Ray_2<Kernel>& r1,
const CGAL::Ray_2<Kernel>& r2);

/*!
returns `true`, if `s1` and `s2` are parallel or if one
of those (or both) is degenerate.
*/
template <typename Kernel>
bool parallel(const CGAL::Segment_2<Kernel>& s1,
const CGAL::Segment_2<Kernel>& s2);

/*!
returns `true`, if `l1` and `l2` are parallel or if one
of those (or both) is degenerate.
*/
template <typename Kernel>
bool parallel(const CGAL::Line_3<Kernel>& l1,
const CGAL::Line_3<Kernel>& l2);

/*!
returns `true`, if `h1` and `h2` are parallel or if one
of those (or both) is degenerate.
*/
template <typename Kernel>
bool parallel(const CGAL::Plane_3<Kernel>& h1,
const CGAL::Plane_3<Kernel>& h2);

/*!
returns `true`, if `r1` and `r2` are parallel or if one
of those (or both) is degenerate.
*/
template <typename Kernel>
bool parallel(const CGAL::Ray_3<Kernel>& r1,
const CGAL::Ray_3<Kernel>& r2);

/*!
returns `true`, if `s1` and `s2` are parallel or if one
of those (or both) is degenerate.
*/
template <typename Kernel>
bool parallel(const CGAL::Segment_3<Kernel>& s1,
const CGAL::Segment_3<Kernel>& s2);

/// @}

/// \defgroup radical_plane_grp CGAL::radical_plane()
/// \ingroup kernel_global_function
/// @{

/*!  returns the radical plane of the two spheres.
 \pre s1 and s2 are not cocentric.
*/
CGAL::Plane_3<Kernel> radical_plane(const CGAL::Sphere_3<Kernel>& s1,
                                    const CGAL::Sphere_3<Kernel>& s2);

/// @}

/// \defgroup radical_line_grp CGAL::radical_line()
/// \ingroup kernel_global_function

/// @{

/*!
returns the radical line of the two circles. 
\pre `c1` and `c2` are not cocentric.
*/
template <typename Kernel>
CGAL::Line_2<Kernel> radical_line(const CGAL::Circle_2<Kernel>& c1,
const CGAL::Circle_2<Kernel>& c2);

/// @}

// Same reason as do_intersect.

/// \defgroup rational_rotation_approximation_grp CGAL::rational_rotation_approximation()
/// \ingroup kernel_global_function

/// \defgroup right_turn_grp CGAL::right_turn()
/// \ingroup kernel_global_function
/// \sa `collinear_grp`
/// \sa `left_turn_grp`
/// \sa `orientation_grp`

/// @{

/*!
returns `true` iff `p`, `q`, and `r` form a right turn.
*/
template <typename Kernel>
bool right_turn(const CGAL::Point_2<Kernel> &p,
const CGAL::Point_2<Kernel> &q,
const CGAL::Point_2<Kernel> &r);

/// @}



/// \defgroup scalar_product_grp CGAL::scalar_product()
/// \ingroup kernel_global_function
/// @{

/*!
returns the scalar product of `u` and `v`.
*/
template <typename Kernel>
Kernel::FT scalar_product( const CGAL::Vector_2<Kernel>& u,
                           const CGAL::Vector_2<Kernel>& v );

/*!
returns the scalar product of `u` and `v`.
*/
template <typename Kernel>
Kernel::FT scalar_product( const CGAL::Vector_3<Kernel>& u,
                           const CGAL::Vector_3<Kernel>& v );
/// @}

/// \defgroup side_of_bounded_circle_grp CGAL::side_of_bounded_circle()
/// \ingroup kernel_global_function
/// \sa `coplanar_side_of_bounded_circle_grp`
/// \sa `side_of_oriented_circle_grp`

/// @{

/*!
returns the relative position of point `t`
to the circle defined by `p`, `q` and `r`. The order
of the points `p`, `q` and `r` does not matter.
\pre `p, q` and `r` are not collinear.
*/
template <typename Kernel>
Bounded_side side_of_bounded_circle(
const CGAL::Point_2<Kernel> &p, 
const CGAL::Point_2<Kernel> &q,
const CGAL::Point_2<Kernel> &r, 
const CGAL::Point_2<Kernel> &t);

/*!
returns the position of the point `t` relative to the circle
that has `pq` as its diameter.
*/
template <typename Kernel>
Bounded_side side_of_bounded_circle(
const CGAL::Point_2<Kernel> &p, 
const CGAL::Point_2<Kernel> &q,
const CGAL::Point_2<Kernel> &t);

/// @}



/// \defgroup side_of_bounded_sphere_grp CGAL::side_of_bounded_sphere()
/// \ingroup kernel_global_function
/// \sa `side_of_oriented_sphere_grp`

/// @{

/*!
returns the relative position of point `t`
to the sphere defined by `p`, `q`, `r`, and `s`. The order
of the points `p`, `q`, `r`, and `s` does not matter.
\pre `p, q, r` and `s` are not coplanar.
*/
template <typename Kernel>
Bounded_side side_of_bounded_sphere(
const CGAL::Point_3<Kernel> &p, 
const CGAL::Point_3<Kernel> &q,
const CGAL::Point_3<Kernel> &r, 
const CGAL::Point_3<Kernel> &s, 
const CGAL::Point_3<Kernel> &t);

/*!
returns the position of the point `t` relative to the sphere
passing through `p`, `q`, and `r` and whose center is in the plane defined
by these three points.
*/
template <typename Kernel>
Bounded_side side_of_bounded_sphere(
const CGAL::Point_3<Kernel> &p, 
const CGAL::Point_3<Kernel> &q,
const CGAL::Point_3<Kernel> &r, 
const CGAL::Point_3<Kernel> &t);

/*!
returns the position of the point `t` relative to the sphere
that has `pq` as its diameter.
*/
template <typename Kernel>
Bounded_side side_of_bounded_sphere(
const CGAL::Point_3<Kernel> &p, 
const CGAL::Point_3<Kernel> &q,
const CGAL::Point_3<Kernel> &t);

/// @}

/// \defgroup side_of_oriented_circle_grp CGAL::side_of_oriented_circle()
/// \ingroup kernel_global_function
/// \sa `side_of_bounded_circle_grp`

/// @{

/*!
returns the relative position of point `test`
to the oriented circle defined by `p`, `q` and `r`.
The order of the points `p`, `q` and `r` is important,
since it determines the orientation of the implicitly
constructed circle.

If `p`, `q` and `r` are collinear, the circle degenerates in a line.
`CGAL::ON_ORIENTED_BOUNDARY` is returned if `test` is also collinear or if two
points are identical, 
otherwise, `side_of_oriented_circle(r, q, test, p)` is returned.

*/
template <typename Kernel>
Oriented_side side_of_oriented_circle(
const CGAL::Point_2<Kernel> &p, 
const CGAL::Point_2<Kernel> &q,
const CGAL::Point_2<Kernel> &r, 
const CGAL::Point_2<Kernel> &test);

/// @}



/// \defgroup side_of_oriented_sphere_grp CGAL::side_of_oriented_sphere()
/// \ingroup kernel_global_function
/// \sa `side_of_bounded_sphere_grp`

/// @{

/*!
returns the relative position of point `test` to the oriented sphere defined
by `p`, `q`, `r` and `s`. The order of the points `p`, `q`, `r`, and `s` is important,
since it determines the orientation of the implicitly constructed
sphere. If the points `p`, `q`, `r` and `s` are positive oriented, positive side
is the bounded interior of the sphere.

In case of degeneracies, `CGAL::ON_ORIENTED_BOUNDARY` is returned
if all points are coplanar. Otherwise, there is a cyclic permutation of the five points
that puts four non coplanar points first, it is used to answer the predicate:
e.g. `CGAL::side_of_oriented_sphere(q, r, s, test, p)` is returned if `q`, `r`, `s`,
and `test` are non coplanar. 
*/
template <typename Kernel>
Oriented_side side_of_oriented_sphere(
const CGAL::Point_3<Kernel> &p, 
const CGAL::Point_3<Kernel> &q,
const CGAL::Point_3<Kernel> &r, 
const CGAL::Point_3<Kernel> &s, 
const CGAL::Point_3<Kernel> &test);

/// @}

/// \defgroup squared_area_grp CGAL::squared_area()
/// \ingroup kernel_global_function

/// @{

/*!
returns the squared area of the triangle defined by the points `p`,
`q` and `r`. 
*/
template <typename Kernel>
Kernel::FT squared_area(const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r);

/// @}

// The same reason as do_intersect.

/// \defgroup squared_distance_grp CGAL::squared_distance()
/// \ingroup kernel_global_function

/// \defgroup squared_radius_grp CGAL::squared_radius()
/// \ingroup kernel_global_function
/// \sa `Circle_2<Kernel>_grp`
/// \sa `Circle_3<Kernel>_grp`
/// \sa `Sphere_3<Kernel>_grp`

/// @{

/*!
compute the squared radius of the circle passing through the points
`p`, `q`, and `r`. \pre `p`, `q`, and `r` are not collinear.
*/
template <typename Kernel>
FT
squared_radius( const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q,
const CGAL::Point_2<Kernel>& r);

/*!
compute the squared radius of the smallest circle passing through `p`,
and `q`, i.e.\ one fourth of the squared distance between `p` and `q`.
*/
template <typename Kernel>
FT
squared_radius( const CGAL::Point_2<Kernel>& p,
const CGAL::Point_2<Kernel>& q);

/*!
compute the squared radius of the smallest circle passing through `p`, 
i.e.\ \f$ 0\f$.
*/
template <typename Kernel>
FT
squared_radius( const CGAL::Point_2<Kernel>& p);

/*!
compute the squared radius of the sphere passing through the points `p`,
`q`, `r` and `s`. \pre `p`, `q`, `r` and `s` are not coplanar.
*/
template <typename Kernel>
FT
squared_radius( const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r,
const CGAL::Point_3<Kernel>& s);

/*!
compute the squared radius of the sphere passing through the points `p`,
`q`, and `r` and whose center is in the same plane as those three points.
*/
template <typename Kernel>
FT
squared_radius( const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q,
const CGAL::Point_3<Kernel>& r);

/*!
compute the squared radius of the smallest circle passing through `p`,
and `q`, i.e.\ one fourth of the squared distance between `p` and `q`.
*/
template <typename Kernel>
FT
squared_radius( const CGAL::Point_3<Kernel>& p,
const CGAL::Point_3<Kernel>& q);

/*!
compute the squared radius of the smallest circle passing through `p`, 
i.e.\ \f$ 0\f$.
*/
template <typename Kernel>
FT
squared_radius( const CGAL::Point_3<Kernel>& p);

/// @}

/// \defgroup unit_normal_grp CGAL::unit_normal()
/// \ingroup kernel_global_function
/// @{

/*!
computes the unit normal vector for the vectors `q-p` and `r-p`.
\pre The points `p`, `q`, and `r` must not be collinear.
*/
CGAL::Vector_3<Kernel> unit_normal( const CGAL::Point_3<Kernel>& p, const CGAL::Point_3<Kernel>& q, const CGAL::Point_3<Kernel>& r );

/// @}

/// \defgroup volume_grp CGAL::volume()
/// \ingroup kernel_global_function
/// \sa `Tetrahedron_3<Kernel>_grp`

/// @{

/*!
Computes the signed volume of the tetrahedron defined by the four points
`p0`, `p1`, `p2` and `p3`.
*/

template <typename Kernel>
Kernel::FT volume(const CGAL::Point_3<Kernel> & p0, const CGAL::Point_3<Kernel> & p1,
                   const CGAL::Point_3<Kernel> & p2, const CGAL::Point_3<Kernel> & p3);

/// @}


/// \defgroup x_equal_grp CGAL::x_equal()
/// \ingroup kernel_global_function
/// \sa `compare_x_grp`
/// \sa `y_equal_grp`
/// \sa `z_equal_grp`

/// @{

/*!
returns `true`, iff `p` and `q`
have the same `x`-coordinate.
*/
template <typename Kernel>
bool x_equal(const CGAL::Point_2<Kernel> &p,
const CGAL::Point_2<Kernel> &q);

/*!
returns `true`, iff `p` and `q`
have the same `x`-coordinate.
*/
template <typename Kernel>
bool x_equal(const CGAL::Point_3<Kernel> &p,
const CGAL::Point_3<Kernel> &q);

/// @}



/// \defgroup y_equal_grp CGAL::y_equal()
/// \ingroup kernel_global_function
/// \sa `compare_y_grp`
/// \sa `x_equal_grp`
/// \sa `z_equal_grp`

/// @{

/*!
returns `true`, iff `p` and `q`
have the same `y`-coordinate.
*/
template <typename Kernel>
bool y_equal(const CGAL::Point_2<Kernel> &p,
const CGAL::Point_2<Kernel> &q);

/*!
returns `true`, iff `p` and `q`
have the same `y`-coordinate.
*/
template <typename Kernel>
bool y_equal(const CGAL::Point_3<Kernel> &p,
const CGAL::Point_3<Kernel> &q);

/// @}


/// \defgroup z_equal_grp CGAL::z_equal()
/// \ingroup kernel_global_function
/// \sa `compare_z_grp`
/// \sa `x_equal_grp`
/// \sa `y_equal_grp`

/// @{

/*!
returns `true`, iff `p` and `q`
have the same `z`-coordinate.
*/
template <typename Kernel>
bool z_equal(const CGAL::Point_3<Kernel> &p,
const CGAL::Point_3<Kernel> &q);

/// @}


/// \defgroup Kernel_operator_plus  CGAL::operator+
/// \ingroup kernel_global_function

/// \defgroup Kernel_operator_minus  CGAL::operator-
/// \ingroup kernel_global_function

/// \defgroup Kernel_operator_prod CGAL::operator*
/// \ingroup kernel_global_function

/// \defgroup do_overlap_grp CGAL::do_overlap()
/// \ingroup kernel_global_function

} /* namespace CGAL */
