
/*!
\ingroup PkgPeriodic2Triangulation2Concepts
\cgalConcept

The concept `Periodic_2TriangulationTraits_2` is the first template parameter of the classes 
`Periodic_2_triangulation_2<Traits, Tds>`. This concept provides the types of 
the geometric primitives used in the triangulation and some function 
object types for the required predicates on those primitives. 

It refines the concept 
`TriangulationTraits_2` from the \cgal \ref PkgTriangulation2 package. It redefines the 
geometric objects, predicates and constructions to work with 
point-offset pairs. In most cases the offsets will be (0,0) and the 
predicates from `TriangulationTraits_2` can be used 
directly. For efficiency reasons we maintain for each functor the 
version without offsets. 

\refines ::TriangulationTraits_2 
In addition to the requirements described for the traits class 
::TriangulationTraits_2, the geometric traits class of a 
Periodic triangulation must fulfill the following 
requirements: 

\hasModel CGAL::Periodic_2_triangulation_traits_2 

\sa `TriangulationTraits_2` 
\sa `CGAL::Periodic_2_triangulation_2<Traits,Tds>` 

*/

class Periodic_2TriangulationTraits_2 {
public:

/// \name Types 
/// @{

// TODO(NGHK): Check
/*! 
The point type. It must be a model of 
`Kernel::Point_2`. 
*/ 
typedef Hidden_type Point_2; 

// TODO(NGHK): Check
/*! 
The segment type. It must be a model 
of `Kernel::Segment_2`. 
*/ 
typedef Hidden_type Segment_2; 

// TODO(NGHK): Check
/*! 
The vector type. It must be a model 
of `Kernel::Vector_2`. 
*/ 
typedef Hidden_type Vector_2; 

// TODO(NGHK): Check
/*! 
The triangle type. It must be a 
model of `Kernel::Triangle_2`. 
*/ 
typedef Hidden_type Triangle_2; 

// TODO(NGHK): Check
/*! 
A type representing an 
axis-aligned rectangle. It must be a model of 
`Kernel::Iso_rectangle_2`. 
*/ 
typedef Hidden_type Iso_rectangle_2; 

// TODO(NGHK): Check
/*! 
The offset type. It must 
be a model of the concept `Periodic_2Offset_2`. 
*/ 
typedef Hidden_type Periodic_2_offset_2; 

/// @} 

/// \name Predicate types 
/// @{

// TODO(NGHK): Check
/*! 

A predicate object that must provide the function operators 

`Comparison_result operator()(Point_2 p, Point_2 q)`, 

which returns `EQUAL` if the \f$ x\f$-coordinates of the two points are equal and 

`Comparison_result operator()(Point_2 p, Point_2 q, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q)`, 

which returns `EQUAL` if the \f$ x\f$-coordinates and \f$ x\f$-offsets of 
the two point-offset pairs are equal. Otherwise it must return a 
consistent order for any two points. \pre `p`, `q` lie inside the domain. 
*/ 
typedef Hidden_type Compare_x_2; 

// TODO(NGHK): Check
/*! 

A predicate object that must provide the function operators 

`Comparison_result operator()(Point_2 p, Point_2 q)`, 

which returns `EQUAL` if the \f$ y\f$-coordinates of the two points are equal and 

`Comparison_result operator()(Point_2 p, Point_2 q, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q)`, 

which returns `EQUAL` if the \f$ y\f$-coordinates and \f$ y\f$-offsets of 
the two point-offset pairs are equal. Otherwise it must return a 
consistent order for any two points. \pre `p`, `q` lie inside the domain. 
*/ 
typedef Hidden_type Compare_y_2; 

// TODO(NGHK): Check
/*! 

Predicate object. Provides the operators: 

`bool operator()(Point p, Point q)` and 

`bool operator()(Point p, Point q, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q)` 

which returns `true` if `p` is before `q` 
according to the \f$ x\f$-ordering of points. 

This predicate is only necessary if the insert function with a range 
of points (using Hilbert sorting) is used. 
*/ 
typedef Hidden_type Less_x_2; 

// TODO(NGHK): Check
/*! 
Predicate object. Provides 
the operators: 

`bool operator()(Point p, Point q)` and 

`bool operator()(Point p, Point q, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q)` 

which returns `true` if `p` is before `q` 
according to the \f$ y\f$-ordering of points. 

This predicate is only necessary if the insert function with a range of 
points (using Hilbert sorting) is used. 

*/ 
typedef Hidden_type Less_y_2; 

// TODO(NGHK): Check
/*! 
A predicate object that must provide the function operators 

`Orientation operator()(Point_2 p, Point_2 q, Point_2 r)`, 

which returns `LEFT_TURN`, `RIGHT_TURN` or `COLLINEAR` 
depending on \f$ r\f$ being, with respect to the oriented line `pq`, 
on the left side, on the right side or on the line. 
and 

`Orientation operator()(Point_2 p, Point_2 q, Point_2 r, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q, Periodic_2_offset_2 o_r)`, 

which returns `LEFT_TURN`, `RIGHT_TURN` or `COLLINEAR` 
depending on `(r,o_r)` being, with respect to the oriented line 
defined by `(p,o_p)(q,o_q)` on the left side, on the right side 
or on the line. 
*/ 
typedef Hidden_type Orientation_2; 

/// @} 

/// \name Constructor types: 
/// Note that the traits must provide exact constructions in order to
/// guarantee exactness of the following construction functors.
/// @{

// TODO(NGHK): Check
/*! 
A constructor object for 
`Point_2`. Provides: 

`Point_2 operator()(Point_2 p,Periodic_2_offset_2 p_o)`, 

which constructs a point from a point-offset pair. 
\pre `p` lies inside the domain. 

*/ 
typedef Hidden_type Construct_point_2; 

// TODO(NGHK): Check
/*! 
A constructor object for 
`Segment_2`. Provides: 

`Segment_2 operator()(Point_2 p,Point_2 q)`, 

which constructs a segment from two points and 

`Segment_2 operator()(Point_2 p,Point_2 q, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q)`, 

which constructs a segment from the points `(p,o_p)` and `(q,o_q)`. 

*/ 
typedef Hidden_type Construct_segment_2; 

// TODO(NGHK): Check
/*! 
A constructor object for 
`Triangle_2`. Provides: 

`Triangle_2 operator()(Point_2 p,Point_2 q,Point_2 r )`, 

which constructs a triangle from three points and 

`Triangle_2 operator()(Point_2 p,Point_2 q,Point_2 r, Periodic_2_offset_2 o_p, Periodic_2_offset_2 o_q, Periodic_2_offset_2 o_r)`, 

which constructs a triangle from the three points `(p,o_p)`, 
`(q,o_q)` and `(r,o_r)`. 
*/ 
typedef Hidden_type Construct_triangle_2; 

/// @} 

/// \name Creation 
/// Only a default constructor, copy constructor and an assignment
/// operator are required. Note that further constructors can be
/// provided.
/// @{

// TODO(NGHK): Check
/*! 
default constructor. 
*/ 
TriangulationTraits_2(); 

// TODO(NGHK): Check
/*! 
Copy constructor 
*/ 
TriangulationTraits_2(TriangulationTraits_2 gtr); 

// TODO(NGHK): Check
/*! 
Assignment operator. 
*/ 
TriangulationTraits_2 operator=(TriangulationTraits_2 gtr); 

/// @} 

/// \name Predicate functions 
/// The following functions give access to the predicate and constructor objects.
/// @{

// TODO(NGHK): Check
/*! 

*/ 
Compare_x_2 compare_x_2_object(); 

// TODO(NGHK): Check
/*! 

*/ 
Compare_y_2 compare_y_2_object(); 

// TODO(NGHK): Check
/*! 

*/ 
Less_x_2 less_x_2_object(); 

// TODO(NGHK): Check
/*! 

*/ 
Less_y_2 less_y_2_object(); 

// TODO(NGHK): Check
/*! 

*/ 
Orientation_2 orientation_2_object(); 

// TODO(NGHK): Check
/*! 

*/ 
Construct_point_2 construct_point_2_object(); 

// TODO(NGHK): Check
/*! 

*/ 
Construct_segment_2 construct_segment_2_object(); 

// TODO(NGHK): Check
/*! 

*/ 
Construct_triangle_2 construct_triangle_2_object(); 

/// @} 

/// \name Access Functions 
/// @{

// TODO(NGHK): Check
/*! 
Sets the 
fundamental domain. This is necessary to evaluate predicates 
correctly. 
\pre `domain` represents a square. 
*/ 
void set_domain(Iso_rectangle_2 domain); 

// TODO(NGHK): Check
/*! 
Returns the 
fundamental domain. 
*/ 
Iso_rectangle_2 get_domain() const; 

/// @}

}; /* end Periodic_2TriangulationTraits_2 */

