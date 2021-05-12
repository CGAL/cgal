
/*!
\ingroup PkgPolytopeDistanceDConcepts
\cgalConcept

This concept defines the requirements for traits classes of
`Width_3<Traits>`.

\cgalHeading{Operations}

Whatever the coordinates of the points are, it is required for the
width-algorithm to have access to the homogeneous representation of
points.

\cgalHasModel CGAL::Width_default_traits_3

\sa `CGAL::Width_3<Traits>`

*/

class WidthTraits_3 {
public:

/// \name Types
/// <I>Notes:</I> If you want to compute the width of a <I>polyhedron</I> then you have to make sure that the point type in the traits class and the point type in the polyhedron class are the same! The same holds for `Traits::Plane_3` and `Polyhedron::Plane_3`.
/// @{

/*!
The point type. The (in)equality tests must be
available. Access to the point coordinates is done via the `get_.()`
functions. Constructing a point is done with the `make_point( )`
operation.
*/
typedef unspecified_type Point_3;

/*!
The plane type. Access to the coefficients of the
plane is made via the `get_.()` functions. Constructing a plane is
done with the `make_plane()` operation.
*/
typedef unspecified_type Plane_3;

/*!
The vector type. There is no need to access the
coefficients of a vector; only constructing is required and is done with
the `make_vector` operation.
*/
typedef unspecified_type Vector_3;

/*!
The traits class for using the convex hull
algorithm. It must be a model of the concept ConvexHullTraits_3.
This class is used only if the width is computed from a set
of points.
*/
typedef unspecified_type ChullTraits;

/*!
Ring type numbers. Internally all numbers are treated as
ring type numbers, i.e., neither \f$ /\f$-operator nor \f$ \sqrt{.}\f$ nor other
inexact operations are used. But because the algorithm does not use any
divisions, but multiplication instead, the numbers can get really big.
Therefore it is recommended to use a ring type number, that provides
values of arbitrary length. Furthermore it is assumed that the underlying
number type of `Point_3`, `Plane_3` and `Vector_3` equals
`RT`.
*/
typedef unspecified_type RT;

/// @}

/// \name Creation
/// Only a default constructor is required.
/// @{

/*!

*/
WidthTraits_3( );

/*!
returns the
homogeneous \f$ x\f$-coordinate of point \f$ p\f$.
*/
RT get_hx(const Point_3& p) const;

/*!
returns the
homogeneous \f$ y\f$-coordinate of point \f$ p\f$.
*/
RT get_hy(const Point_3& p) const;

/*!
returns the
homogeneous \f$ z\f$-coordinate of point \f$ p\f$.
*/
RT get_hz(const Point_3& p) const;

/*!
returns the
homogenizing coordinate of point \f$ p\f$.
*/
RT get_hw(const Point_3& p) const;

/*!
returns all homogeneous coordinates
of point \f$ p\f$ at once.
*/
void get_point_coordinates(const Point_3& p, RT& px,
RT& py, RT& pz, RT& ph) const;

/*!
returns the first
coefficient of plane \f$ f\f$.
*/
RT get_a(const Plane_3& f) const;

/*!
returns the second
coefficient of plane \f$ f\f$.
*/
RT get_b(const Plane_3& f) const;

/*!
returns the third
coefficient of plane \f$ f\f$.
*/
RT get_c(const Plane_3& f) const;

/*!
returns the fourth
coefficient of plane \f$ f\f$.
*/
RT get_d(const Plane_3& f) const;

/*!
returns all four plane coefficients of \f$ f\f$ at once.
*/
void get_plane_coefficients(const Plane_3& f,
RT& a, RT& b, RT& c, RT& d)
const;

/// @}

}; /* end WidthTraits_3 */

/*!
returns a point of type
`Point_3` with homogeneous coordinates \f$ hx\f$, \f$ hy\f$, \f$ hz\f$ and \f$ hw\f$.
\relates WidthTraits_3
*/
Point_3 make_point(const RT& hx, const RT& hy, const RT& hz,
const RT& hw) const;

/*!
returns a plane of type `Plane_3`
with coefficients \f$ a\f$, \f$ b\f$, \f$ c\f$ and \f$ d\f$.
\relates WidthTraits_3
*/
Plane_3 make_plane(const RT& a, const RT& b, const
RT& c, const RT& d) const;

/*!
returns a vector of type `Vector_3` with the four
homogeneous coefficients \f$ a\f$, \f$ b\f$, \f$ c\f$ and 1.
\relates WidthTraits_3
*/
Vector_3 make_vector(const RT& a, const RT& b, const
RT& c) const;
