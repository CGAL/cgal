
/*!
\ingroup PkgPeriodic3Triangulation3Concepts
\cgalConcept

The concept `Periodic_3DelaunayTriangulationTraits_3` is the first template parameter of the classes 
`Periodic_3_Delaunay_triangulation_3` and 
`Periodic_3_triangulation_3`. It refines the concept 
`DelaunayTriangulationTraits_3` from the \cgal \ref PkgTriangulation3 package. 
It redefines the geometric objects, predicates and constructions to 
work with point-offset pairs. In most cases the offsets will be 
(0,0,0) and the predicates from `DelaunayTriangulationTraits_3` 
can be used directly. For efficiency reasons we maintain for each 
functor the version without offsets. 

\cgalRefines `DelaunayTriangulationTraits_3` 

\cgalHasModel CGAL::Periodic_3_triangulation_traits_3 

\sa `DelaunayTriangulationTraits_3` 

In addition to the requirements described for the traits class 
DelaunayTriangulationTraits_3, the geometric traits class of a 
Periodic Delaunay triangulation must fulfill the following 
requirements.

\note The optional types must be provided in any case, however they
can be replaced by dummy types if the respective functions are not
used.
*/
class Periodic_3DelaunayTriangulationTraits_3 {
public:

/// \name Types 
/// @{

/*! 
The point type. It must be a model of `Kernel::Point_3`. 
*/ 
typedef Hidden_type Point_3; 

/*! 
The vector type. It must be a model of 
`Kernel::Vector_3`. 
*/ 
typedef Hidden_type Vector_3; 

/*! 
The offset type. It must be a 
model of the concept `Periodic_3Offset_3`. 
*/ 
typedef Hidden_type Periodic_3_offset_3; 

/*! 
A type representing an axis-aligned 
cuboid. It must be a model of `Kernel::Iso_cuboid_3`. 
*/ 
typedef Hidden_type Iso_cuboid_3; 

/// @}

/// \name
/// The following three types represent geometric primitives in \f$
/// \mathbb R^3\f$. They are required to provide functions converting
/// primitives from \f$ \mathbb T_c^3\f$ to \f$ \mathbb R^3\f$,
/// i.e.\ constructing representatives in \f$ \mathbb R^3\f$.
/// @{

/*! 
A segment type. It must be a model of `Kernel::Segment_3`. 
*/ 
typedef Hidden_type Segment_3; 

/*! 
A triangle type. It must be a model of 
`Kernel::Triangle_3`. 
*/ 
typedef Hidden_type Triangle_3; 

/*! 
A tetrahedron type. It must be a model 
of `Kernel::Tetrahedron_3`. 
*/ 
typedef Hidden_type Tetrahedron_3; 

/*! 
A predicate object that must provide the function operators 

`Comparison_result operator()(Point_3 p, Point_3 q)`, 

which returns `EQUAL` if the two points are equal and 

`Comparison_result operator()(Point_3 p, Point_3 q, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q)`, 

which returns `EQUAL` if the two point-offset pairs are equal. 
Otherwise it must return a consistent order for any two points chosen 
in a same line. 
\pre `p`, `q` lie inside the domain. 
*/ 
typedef Hidden_type Compare_xyz_3; 

/*! 
A predicate object that must provide the function operators 

`Orientation operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`, 

which returns `POSITIVE`, if `s` lies on the positive side of 
the oriented plane `h` defined by `p`, `q`, and `r`, 
returns `NEGATIVE` if `s` lies on the negative side of 
`h`, and returns `COPLANAR` if `s` lies on `h` and 

`Orientation operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_s)`, 

which returns `POSITIVE`, if the point-offset pair `(s,o_s)` 
lies on the positive side of the oriented plane `h` defined by 
`(p,o_p)`, `(q,o_q)`, and `(r,o_r)`, 
returns `NEGATIVE` if `(s,o_s)` lies on the negative side of 
`h`, and returns `COPLANAR` if `(s,o_s)` lies on `h`. 
\pre `p`, `q`, `r`, `s` lie inside the domain. 
*/ 
typedef Hidden_type Orientation_3; 

/*! 
A predicate object that must provide the function operators 

`Oriented_side operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Point_3 t)`, 

which determines on which side of the oriented sphere circumscribing 
`p, q, r, s` the point `t` lies and 

`Oriented_side operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Point_3 t, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_s, Periodic_3_offset_3 o_t)`, 

which determines on which side of the oriented sphere circumscribing 
`(p,o_p), (q,o_q), (r,o_r), (s,o_s)` the point-offset pair 
`(t,o_t)` lies. 
\pre `p`, `q`, `r`, `s`, `t` lie inside the domain. 
*/ 
typedef Hidden_type Side_of_oriented_sphere_3; 

/*! 
A predicate object that must provide the function operators 

`Comparison_result operator()(Point_3 p, Point_3 q, Point_3 r)`, 

which compares the distance between `p` and `q` to the distance 
between `p` and `r` and 

`Comparison_result operator()(Point_3 p, Point_3 q, Point_3 r, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r)`, 

which compares the distance between `(p,o_p)` and `(q,o_q)` to 
the distance between `(p,o_p)` and `(r,o_r)`. 
\pre `p`, `q`, `r` lie inside the domain. 
*/ 
typedef Hidden_type Compare_distance_3; 

/// @}

/// \name
/// In addition, only when vertex removal is used, the traits class must provide the following predicate objects
/// @{

/*! 
A predicate object that must provide the function operators 

`Orientation operator()(Point_3 p, Point_3 q, Point_3 r)`, 

which returns `COLLINEAR`, if the points are collinear; otherwise 
it must return a consistent orientation for any three points chosen in 
a same plane and 

`Orientation operator()(Point_3 p, Point_3 q, Point_3 r Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r)`, 

which returns `COLLINEAR`, if the point-offset pairs are 
collinear; otherwise it must return a consistent orientation for any 
three point-offset pairs chosen in a same plane. 
\pre `p`, `q`, `r` lie inside the domain. 
*/ 
typedef Hidden_type Coplanar_orientation_3; 

/*! 
A predicate object that must provide the function operators 

`Bounded_side operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`, 

which determines the bounded side of the circle defined by `p, q`, 
and `r` on which the point `s` lies and 

`Bounded_side operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_s)`, 

which determines the bounded side of the circle defined by 
`(p,o_p), (q,o_q)`, and `(r,o_r)` on which the point-offset pair 
`(s,o_s)` lies. 
\pre `p,q,r`, and `s` are coplanar and `p,q`, and `r` are not collinear, `(p,o_p),(q,o_q),(r,o_r)`, and `(s,o_s)` are coplanar and `(p,o_p),(q,o_q)`, and `(r,o_r)` are not collinear, respectively, and `p`, `q`, `r`, `s`, `t` lie inside the domain. 
*/ 
typedef Hidden_type Coplanar_side_of_bounded_circle_3; 

/// @}

/// \name
/// In addition, only when `is_Gabriel` is used, the traits class must
/// provide the following predicate object:
/// @{

/*! 
A predicate object that must provide the function operators 

`Bounded_side operator()(Point_3 p, Point_3 q, Point_3 t)`, 

which returns the position of the point `t` relative to the sphere 
that has `pq` as its diameter, 

`Bounded_side operator()(Point_3 p, Point_3 q, Point_3 t, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_t)`, 

which returns the position of the point-offset pair `(t,o_t)` 
relative to the sphere that has `(p,o_p)(q,o_q)` as its diameter, 

`Bounded_side operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 t)`, 

which returns the position of the point `t` relative to the sphere 
passing through `p, q`, and `r` and whose center is in the 
plane defined by these three points, 

`Bounded_side operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 t, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_q)`, 

which returns the position of the point-offset pair `(t,o_t)` 
relative to the sphere passing through `(p,o_p), (q,o_q)`, and 
`(r,o_r)` and whose center is in the plane defined by these three 
point-offset pairs, 

`Bounded_side operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Point_3 t)`, 

which returns the relative position of point `t` to the sphere 
defined by `p, q, r`, and `s`; the order of the points `p, q, r`, and `s` does not matter, and 

`Bounded_side operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Point_3 t, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_s, Periodic_3_offset_3 o_q)`, 

which returns the relative position of the point-offset pair 
`(t,o_t)` to the sphere defined by `(p,o_p), (q,o_q), (r,o_r)`, and `(s,o_s)`; the order of the point-offset pairs 
`(p,o_p), (q,o_q), (r,o_r)`, and `(s,o_s)` does not matter. 
\pre `p, q, r`, and `s` are not coplanar, `(p,o_p), (q,o_q), (r,o_r)`, and `(s,o_s)` are not coplanar, `p`, `q`, `r`, `s`, `t` lie inside the domain. 
*/ 
typedef Hidden_type Side_of_bounded_sphere_3; 

/// @}

/// \name 
/// Note that the traits must provide exact constructions in order to
/// guarantee exactness of the following construction functors.
/// @{

/*! 
A constructor object that must provide the function operator 

`Point_3 operator()(Point_3 p, Periodic_3_offset_3 o_p)`, 

which constructs a point from a point-offset pair. 
\pre `p` lies inside the domain. 
*/ 
typedef Hidden_type Construct_point_3; 

/*! 
A constructor object that must provide the function operators 

`Segment_3 operator()(Point_3 p, Point_3 q)`, 

which constructs a segment from two points and 

`Segment_3 operator()(Point_3 p, Point_3 q, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q)`, 

which constructs a segment from two point-offset pairs. 
\pre `p`, `q` lie inside the domain. 
*/ 
typedef Hidden_type Construct_segment_3; 

/*! 
A constructor object that must provide the function operators 

`Triangle_3 operator()(Point_3 p, Point_3 q, Point_3 r )`, 

which constructs a triangle from three points and 

`Triangle_3 operator()(Point_3 p, Point_3 q, Point_3 r, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r)`, 

which constructs a triangle from three point-offset pairs. 
\pre `p`, `q`, `r` lie inside the domain. 
*/ 
typedef Hidden_type Construct_triangle_3; 

/*! 
A constructor object that must provide the function operators 

`Tetrahedron_3 operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`, 

which constructs a tetrahedron from four points and 

`Tetrahedron_3 operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_s)`, 

which constructs a tetrahedron from four point-offset pairs. 
\pre `p`, `q`, `r`, `s` lie inside the domain. 
*/ 
typedef Hidden_type Construct_tetrahedron_3; 

/// @}

/// \name
/// In addition, only when the dual operations are used, the traits
/// class must provide the following constructor object:
/// @{

/*! 
A constructor object that must provide the function operators 

`Point_3 operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`, 

which constructs the circumcenter of four points and 

`Point_3 operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_s)`, 

which constructs the circumcenter of four point-offset pairs. 
\pre `p`, `q`, `r` and `s` as well as `(p,o_p)`, `(q,o_q)`, `(r,o_r)` and `(s,o_s)` must be non coplanar. `p`, `q`, `r`, `s` lie inside the domain. 
*/ 
typedef Hidden_type Construct_circumcenter_3; 

/// @} 

/// \name Creation 
/// @{

/*! 
Default constructor. 
*/ 
Periodic_3_triangulation_traits_3(); 

/*! 
Copy constructor. 
*/ 
Periodic_3_triangulation_traits_3(const Periodic_triangulation_traits_3 & tr); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
Set the size of the 
fundamental domain. This is necessary to evaluate predicates 
correctly. 
\pre `domain` represents a cube. 
*/ 
void set_domain(Iso_cuboid_3 domain); 

/// @} 

/// \name Operations 
/// The following functions give access to the predicate and construction objects:
/// @{

/*! 

*/ 
Compare_xyz_3 compare_xyz_3_object(); 

/*! 

*/ 
Orientation_3 orientation_3_object(); 

/*! 

*/ 
Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object(); 

/*! 

*/ 
Compare_distance_3 compare_distance_3_object(); 

/// @}

/// \name
/// The following functions must be provided if vertex removal is
/// used; otherwise dummy functions can be provided.
/// @{

/*! 

*/ 
Coplanar_orientation_3 coplanar_3_orientation_3_object(); 

/*! 

*/ 
Coplanar_side_of_bounded_circle_3 
coplanar_side_of_bounded_circle_3_object(); 

/// @}

/// \name
/// The following function must be provided only if the `is_Gabriel`
/// methods of `Periodic_3_Delaunay_triangulation_3` are used;
/// otherwise a dummy function can be provided.
/// @{

/*! 

*/ 
Side_of_bounded_sphere_3 side_of_bounded_sphere_3_object(); 

/*! 

*/ 
Construct_segment_3 construct_segment_3_object(); 

/*! 

*/ 
Construct_triangle_3 construct_triangle_3_object(); 

/*! 

*/ 
Construct_tetrahedron_3 construct_tetrahedron_3_object(); 

/// @}

/// \name 
/// The following function must be provided only if the methods of
/// `Periodic_3_Delaunay_triangulation_3` returning elements of the
/// Voronoi diagram are used; otherwise a dummy function can be
/// provided:
/// @{

/*! 

*/ 
Construct_circumcenter_3 construct_circumcenter_3_object(); 

/// @}

}; /* end Periodic_3DelaunayTriangulationTraits_3 */

