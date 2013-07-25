
/*!
\ingroup PkgNef2Concepts
\cgalConcept

`ExtendedKernelTraits_2` is a kernel concept providing extended
geometry\cgalFootnote{It is called extended geometry for simplicity,
though it is not a real geometry in the classical sense}. Let `K` be
an instance of the data type `ExtendedKernelTraits_2`. The central
notion of extended geometry are extended points. An extended point
represents either a standard affine point of the Cartesian plane or a
non-standard point representing the equivalence class of rays where
two rays are equivalent if one is contained in the other.

Let \f$ R\f$ be an infinimaximal number\cgalFootnote{A finite but very large number.}, \f$ F\f$ be the square box with corners \f$ NW(-R,R)\f$, \f$ NE(R,R)\f$, 
\f$ SE(R,-R)\f$, and \f$ SW(-R,-R)\f$. Let \f$ p\f$ be a non-standard point and let 
\f$ r\f$ be a ray defining it. If the frame \f$ F\f$ contains the source point 
of \f$ r\f$ then let \f$ p(R)\f$ be the intersection of \f$ r\f$ with the frame \f$ F\f$, 
if \f$ F\f$ does not contain the source of \f$ r\f$ then \f$ p(R)\f$ is undefined. 
For a standard point let \f$ p(R)\f$ be equal to \f$ p\f$ if \f$ p\f$ is contained in 
the frame \f$ F\f$ and let \f$ p(R)\f$ be undefined otherwise. Clearly, for any 
standard or non-standard point \f$ p\f$, \f$ p(R)\f$ is defined for any 
sufficiently large \f$ R\f$. Let \f$ f\f$ be any function on standard points, 
say with \f$ k\f$ arguments. We call \f$ f\f$ <I>extensible</I> if for any \f$ k\f$ 
points \f$ p_1\f$, \f$ \ldots\f$ , \f$ p_k\f$ the function value 
\f$ f(p_1(R),\ldots,p_k(R))\f$ is constant for all sufficiently large 
\f$ R\f$. We define this value as \f$ f(p_1,\ldots,p_k)\f$. Predicates like 
lexicographic order of points, orientation, and incircle tests are 
extensible. 

An extended segment is defined by two extended points such that it is 
either an affine segment, an affine ray, an affine line, or a segment 
that is part of the square box. Extended directions extend the affine 
notion of direction to extended objects. 

This extended geometry concept serves two purposes. It offers 
functionality for changing between standard affine and extended 
geometry. At the same time it provides extensible geometric primitives 
on the extended geometric objects. 

\cgalHasModel `CGAL::Extended_cartesian<FT>`
\cgalHasModel `CGAL::Extended_homogeneous<RT>`
\cgalHasModel `CGAL::Filtered_extended_homogeneous<RT>`

*/

class ExtendedKernelTraits_2 {
public:

/// \name Affine kernel types 
/// @{

/*!
the standard affine kernel. 

*/ 
typedef unspecified_type Standard_kernel; 

/*!
the standard ring type. 

*/ 
typedef unspecified_type Standard_RT; 

/*!
standard points. 

*/ 
typedef unspecified_type Standard_point_2; 

/*!
standard segments. 

*/ 
typedef unspecified_type Standard_segment_2; 

/*!
standard oriented lines. 

*/ 
typedef unspecified_type Standard_line_2; 

/*!
standard directions. 

*/ 
typedef unspecified_type Standard_direction_2; 

/*!
standard rays. 

*/ 
typedef unspecified_type Standard_ray_2; 

/*!
standard affine transformations. 

*/ 
typedef unspecified_type Standard_aff_transformation_2; 

/// @} 

/// \name Extended kernel types 
/// @{

/*!
the ring type of our extended kernel. 

*/ 
typedef unspecified_type RT; 

/*!
extended points. 

*/ 
typedef unspecified_type Point_2; 

/*!
extended segments. 

*/ 
typedef unspecified_type Segment_2; 

/*!
extended directions. 

*/ 
typedef unspecified_type Direction_2; 

/*!
a type descriptor for extended points. 

*/ 
enum Point_type { SWCORNER, LEFTFRAME, NWCORNER, BOTTOMFRAME, STANDARD, TOPFRAME, SECORNER, RIGHTFRAME, NECORNER }; 

/// @} 

/// \name Interfacing the affine kernel types 
/// @{

/*!

creates an extended point and initializes it to the standard point 
`p`. 
*/ 
Point_2 construct_point(const Standard_point_2& p) ; 

/*!

creates an extended point and initializes it to the equivalence class 
of all the rays underlying the oriented line `l`. 
*/ 
Point_2 construct_point(const Standard_line_2& l) ; 

/*!
creates an extended point and initializes it 
to the equivalence class of all the rays underlying the oriented line 
`l(p1,p2)`. 
*/ 
Point_2 construct_point(const Standard_point_2& p1, const 
Standard_point_2& p2) ; 

/*!
creates an extended point and initializes 
it to the equivalence class of all the rays underlying the ray 
starting in `p` in direction `d`. 
*/ 
Point_2 construct_point(const Standard_point_2& p, const 
Standard_direction_2& d) ; 

/*!
creates an extended point and initializes it to the equivalence 
class of all the rays underlying the oriented line opposite to 
`l`. 
*/ 
Point_2 construct_opposite_point(const Standard_line_2& l) 
; 

/*!
determines the type of 
`p` and returns it. 
*/ 
Point_type type(const Point_2& p) ; 

/*!
returns `true` iff 
`p` is a standard point. 
*/ 
bool is_standard(const Point_2& p) ; 

/*!
returns 
the standard point represented by `p`. \pre `K.is_standard(p)`. 
*/ 
Standard_point_2 standard_point(const Point_2& p) ; 

/*!
returns 
the oriented line representing the bundle of rays defining `p`. 
\pre `!K.is_standard(p)`. 
*/ 
Standard_line_2 standard_line(const Point_2& p) ; 

/*!
a ray 
defining `p`. \pre `!K.is_standard(p)`. 
*/ 
Standard_ray_2 standard_ray(const Point_2& p) ; 

/*!
returns the point on the northeast frame 
corner. 
*/ 
Point_2 NE() ; 

/*!
returns the point on the southeast frame 
corner. 
*/ 
Point_2 SE() ; 

/*!
returns the point on the northwest frame 
corner. 
*/ 
Point_2 NW() ; 

/*!
returns the point on the southwest frame 
corner. 
*/ 
Point_2 SW() ; 

/// @} 

/// \name Geometric kernel calls 
/// @{

/*!
returns the source 
point of `s`. 
*/ 
Point_2 source(const Segment_2& s) ; 

/*!
returns the target 
point of `s`. 
*/ 
Point_2 target(const Segment_2& s) ; 

/*!
constructs a segment `pq`. 
*/ 
Segment_2 construct_segment(const Point_2& p, const Point_2& 
q) ; 

/*!
returns the orientation of `p` with respect to the line through 
`s`. 
*/ 
int orientation(const Segment_2& s, const Point_2& p) 
; 

/*!
returns the orientation of `p3` with respect to 
the line through `p1p2`. 
*/ 
int orientation(const Point_2& p1, const Point_2& p2, const 
Point_2& p3) ; 

/*!
return true iff the `p3` is left of the line 
through `p1p2`. 
*/ 
bool left_turn(const Point_2& p1, const Point_2& p2, const 
Point_2& p3) ; 

/*!
return true iff 
`s` is degenerate. 
*/ 
bool is_degenerate(const Segment_2& s) ; 

/*!
returns the lexicographic order of `p1` and `p2`. 
*/ 
int compare_xy(const Point_2& p1, const Point_2& p2) 
; 

/*!
returns the order on the \f$ x\f$-coordinates of `p1` and `p2`. 

*/ 
int compare_x(const Point_2& p1, const Point_2& p2) 
; 

/*!
returns the order on the \f$ y\f$-coordinates of `p1` and `p2`. 

*/ 
int compare_y(const Point_2& p1, const Point_2& p2) 
; 

/*!
returns the point of intersection of the lines supported by 
`s1` and `s2`. \pre The intersection point exists. 
*/ 
Point_2 intersection( const Segment_2& s1, const Segment_2& 
s2) ; 

/*!
returns the direction of the vector `p2` - 
`p1`. 
*/ 
Direction_2 construct_direction( const Point_2& p1, const 
Point_2& p2) ; 

/*!
returns `true` iff 
`d2` is in the interior of the counterclockwise angular sector 
between `d1` and `d3`. 
*/ 
bool strictly_ordered_ccw(const Direction_2& d1, const 
Direction_2& d2, const Direction_2& d3) ; 

/*!
returns `true` iff `p2` is 
in the relative interior of the segment `p1p3`. 
*/ 
bool strictly_ordered_along_line( const Point_2& p1, const 
Point_2& p2, const Point_2& p3) ; 

/*!
returns true iff `s` contains `p`. 
*/ 
bool contains(const Segment_2& s, const Point_2& p) 
; 

/*!
returns true iff 
\f$ \|p1-p2\| < \|p3-p4\|\f$. 
*/ 
bool first_pair_closer_than_second( const Point_2& p1, const 
Point_2& p2, const Point_2& p3, const Point_2& p4) ; 

/*!
returns a unique 
identifier for kernel object Input/Output. Usually this should be the 
name of the model. 
*/ 
const char* output_identifier() ; 

/// @}

}; /* end ExtendedKernelTraits_2 */

