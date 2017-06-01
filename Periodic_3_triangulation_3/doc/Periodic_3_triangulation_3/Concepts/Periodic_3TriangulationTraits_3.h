
/*!
\ingroup PkgPeriodic3Triangulation3Concepts
\cgalConcept

The concept `Periodic_3TriangulationTraits_3` is the first template parameter of the class
`Periodic_3_triangulation_3`. It refines the concept
`TriangulationTraits_3` from the \cgal 3D Triangulations.
It redefines the geometric objects, predicates and constructions to
work with point-offset pairs. In most cases the offsets will be
(0,0,0) and the predicates from `TriangulationTraits_3`
can be used directly. For efficiency reasons we maintain for each
functor the version without offsets.

\cgalRefines TriangulationTraits_3

\cgalHasModel CGAL::Periodic_3_triangulation_traits_3

\sa Periodic_3DelaunayTriangulationTraits_3
\sa Periodic_3RegularTriangulationTraits_3

In addition to the requirements described for the traits class
TriangulationTraits_3, the geometric traits class of a
Periodic triangulation must fulfill the following
requirements.
*/
class Periodic_3TriangulationTraits_3 {
public:

/// \name Types
/// @{

/*!
The point type. It must be a model of `Kernel::Point_3`.
*/
typedef unspecified_type Point_3;

/*!
The vector type. It must be a model of
`Kernel::Vector_3`.
*/
typedef unspecified_type Vector_3;

/*!
The offset type. It must be a
model of the concept `Periodic_3Offset_3`.
*/
typedef unspecified_type Periodic_3_offset_3;

/*!
A type representing an axis-aligned
cuboid. It must be a model of `Kernel::Iso_cuboid_3`.
*/
typedef unspecified_type Iso_cuboid_3;

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
typedef unspecified_type Segment_3;

/*!
A triangle type. It must be a model of
`Kernel::Triangle_3`.
*/
typedef unspecified_type Triangle_3;

/*!
A tetrahedron type. It must be a model
of `Kernel::Tetrahedron_3`.
*/
typedef unspecified_type Tetrahedron_3;

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
typedef unspecified_type Compare_xyz_3;

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
typedef unspecified_type Orientation_3;

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
typedef unspecified_type Construct_point_3;

/*!
A constructor object that must provide the function operators

`Segment_3 operator()(Point_3 p, Point_3 q)`,

which constructs a segment from two points and

`Segment_3 operator()(Point_3 p, Point_3 q, Periodic_3_offset_3 o_p, Periodic_3_offset_3 o_q)`,

which constructs a segment from two point-offset pairs.
\pre `p`, `q` lie inside the domain.
*/
typedef unspecified_type Construct_segment_3;

/*!
A constructor object that must provide the function operators

`Triangle_3 operator()(Point_3 p, Point_3 q, Point_3 r )`,

which constructs a triangle from three points and

`Triangle_3 operator()(Point_3 p, Point_3 q, Point_3 r, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r)`,

which constructs a triangle from three point-offset pairs.
\pre `p`, `q`, `r` lie inside the domain.
*/
typedef unspecified_type Construct_triangle_3;

/*!
A constructor object that must provide the function operators

`Tetrahedron_3 operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`,

which constructs a tetrahedron from four points and

`Tetrahedron_3 operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_q, Periodic_3_offset_3 o_r, Periodic_3_offset_3 o_s)`,

which constructs a tetrahedron from four point-offset pairs.
\pre `p`, `q`, `r`, `s` lie inside the domain.
*/
typedef unspecified_type Construct_tetrahedron_3;

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
Set the fundamental domain. This is necessary to evaluate predicates correctly.
\pre `domain` represents a cube.
*/
void set_domain(const Iso_cuboid_3& domain);

/*!
Returns the fundamental domain.
*/
const Iso_cuboid_3& get_domain() const;

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

/// @}

/// @{

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

}; /* end Periodic_3TriangulationTraits_3 */

