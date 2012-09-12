/*!


*/
class Arr_algebraic_segment_traits_2::Construct_curve_2 {
public:

/// \name Object Creation Functors 
/// @{

/*! 
Returns a `Curve_2` object that represents the curve defined by 
the polynomial `p` 
*/ 
Curve_2 operator() (Polynomial_2 p); 

/*! 
Returns a `Curve_2` object specified by `s`. 
The passed string represents the defining polynomial of the curve 
and must be given in a MAPLE-readable format using "x" as first 
and "y" as second variable, e.g., 
"(x&circ;3*y-2*x)*(-6*x-y&circ;3*x&circ;6)" 
for integer coefficients, and "3/2*x*y&circ;4-5/7*x&circ;2+3/1" 
for rational coefficients. 
*/ 
Curve_2 operator() (std::string s); 

/// @}

}; /* end Arr_algebraic_segment_traits_2::Construct_curve_2 */

/*!


*/
class Arr_algebraic_segment_traits_2::Construct_point_2 {
public:

/// \name Object Creation Functors 
/// @{

/*! 
Returns a `Point_2` object that represents the `arcno`-th 
point in the fiber of `cv` at \f$ x\f$-coordinate `x`, 
counted from the bottom, starting with zero. 
\pre (`cv` must not have a vertical line at `x`, 
and \f$ 0\leq arcno < c\f$, where \f$ c\f$ is the number of points 
in the fiber of `cv` at `x`.) 
*/ 
Point_2 operator() (Algebraic_real_1 x, Curve_2 cv, int arcno); 

/*! 
Returns a `Point_2` object that represents the 
point on `xcv` at \f$ x\f$-coordinate `x` 
\pre (`x` is in the \f$ x\f$-range of `xcv`.) 
*/ 
Point_2 operator() (Algebraic_real_1 x, X_monotone_curve_2 xcv); 

/*! 
Returns a `Point_2` object that represents (x,y) 
*/ 
Point_2 operator() (Algebraic_real_1 x, Algebraic_real_1 y); 

/*! 
Returns a `Point_2` object that represents (x,y) 
*/ 
Point_2 operator() (Coefficient x, Coefficient y); 

/*! 
Returns a `Point_2` object that represents (x,y) 
*/ 
Point_2 operator() (Bound x, Bound y); 

/*! 
Returns a `Point_2` object that represents (x,y) 
*/ 
Point_2 operator() (int x, int y); 

/// @}

}; /* end Arr_algebraic_segment_traits_2::Construct_point_2 */

/*!


*/
class Arr_algebraic_segment_traits_2::Construct_x_monotone_segment_2 {
public:

/// \name Object Creation Functors 
/// @{

/*! 
Writes a sequence of `X_monotone_curve_2` objects (terminal 
segments) into `out`. These terminal segments compose an 
\f$ x\f$-monotone (or vertical) segment of the curve `cv` that 
starts in `end_min`, and ends in `end_max`. 
\pre (`end_min` must have a unique \f$ x\f$-monotone segment 
to its right, or `end_max` must have a unique 
\f$ x\f$-monotone segment to its left. Furthermore, 
`end_min` and `end_max` must be connected 
by an \f$ x\f$-monotone segment of `cv`) 
*/ 
template<class OutputIterator> OutputIterator 
operator() (Curve_2 cv, Point_2 end_min, Point_2 end_max, 
OutputIterator out); 

/*! 
Writes a sequence of `X_monotone_curve_2` objects into `out`. 
These segments form an \f$ x\f$-monotone (or vertical) 
segment of the curve `cv`. 

If `site_of_p==POINT_IN_INTERIOR`, the maximal segment is 
returned that contains `p` in its interior. 

returned that contains `p` as its left endpoint. 

returned that contains `p` as its left endpoint. 
\pre (If `site_of_p==POINT_IN_INTERIOR`, `p` 
must be an interior point of an \f$ x\f$-monotone or a vertical 
segment. 
must either have a unique \f$ x\f$-monotone segment to the right, 
or a vertical segment from `p` upwards. 
must either have a unique \f$ x\f$-monotone segment to the left, 
or a vertical segment from `p` downwards.) 
*/ 
template<class OutputIterator> OutputIterator 
operator() (Curve_2 cv, Point_2 p, 
Site_of_point site_of_p, 
OutputIterator out); 

/*! 
Writes a sequence of `X_monotone_curve_2` objects into `out`. 
These segments form a straight-line segment connecting 
the points `p` and `q`. If `p` and `q` share the 
same \f$ x\f$-coordinate, the constructed vertical segment consists of 
only one `X_monotone_curve_2` object and can be computed 
efficiently. In the non-vertical case, 
the construction is only possible if `p` and `q` 
have both rational x- and y-coordinates. 
\pre (`p` must not be equal to `q`.) 

*/ 
template<class OutputIterator> OutputIterator 
operator() (Point_2 p, Point_2 q, 
OutputIterator out); 

/// @}

}; /* end Arr_algebraic_segment_traits_2::Construct_x_monotone_segment_2 */

/*!


Models the `ArrangementTraits_2::Curve_2` concept. 
Represents algebraic curves. Internally, the type stores 
topological-geometric information about the particular curve. 
In order to use internal caching, instances should only be created 
using the `Construct_curve_2` functor of the traits class. 

*/
class Arr_algebraic_segment_traits_2::Curve_2 {
public:

/// \name Modifiers 
/// @{

/*! 
returns the defining polynomial of the curve. 
*/ 
Polynomial_2 polynomial () const; 

/// @}

}; /* end Arr_algebraic_segment_traits_2::Curve_2 */

/*!


Models the `ArrangementBasicTraits_2::Point_2` concept. 
Represents points in \f$ \mathbb{R}^2\f$. Intersection points of algebraic curves 
are in general non-rational, so we need a data structure that is capable 
of representing arbitrary points with algebraic coordinates. 

The traits class represents algebraic coordinates by the type 
`Algebraic_real_1`, which is a model of the `AlgebraicReal_1` 
concept. 
A point \f$ p\f$ is stored by a triple \f$ (x,cv,arcno)\f$, 
where \f$ x\f$ is the \f$ x\f$-coordinate of a point, \f$ cv\f$ is an instance 
of `Curve_2` that contains the point, (and has no vertical line at \f$ x\f$), 
and \f$ arcno\f$ is an `int`, denoting that \f$ p\f$ is met as the \f$ arcno\f$-th 
point when shooting a vertical ray at \f$ x\f$, starting from \f$ -\infty\f$ 
(where counting starts with \f$ 0\f$). 

In addition to the methods listed below, the copy constructor and assignment 
operator for `Point_2` objects are also supported. 

The functor `Construct_point_2` constructs `Point_2` instances. 

*/
class Arr_algebraic_segment_traits_2::Point_2 {
public:

/// \name Modifiers 
/// @{

/*! 
returns the \f$ x\f$-coordinate of `p`. 
*/ 
Algebraic_real_1 x () const; 

/*! 
returns the \f$ y\f$-coordinates of `p`. 

<B>Attention:</B> As described above, points are not stored 
by their \f$ y\f$-coordinate in `Algebraic_real_1` representation. In fact, 
this representation must be computed on demand, and might become quite 
costly for points defined by high-degree polynomials. Therefore, it is 
recommended to avoid to call this function as much as possible. 
*/ 
Algebraic_real_1 y () const; 

/*! 
returns a `Curve_2` instance that `p`is part of. 
*/ 
Curve_2 curve () const; 

/*! 
returns the arc number of `p`. 
*/ 
int arcno () const; 

/*! 
returns double-approximations of the \f$ x\f$- and \f$ y\f$-coordinates. 
*/ 
std::pair<double,double> to_double () const; 

/// @}

}; /* end Arr_algebraic_segment_traits_2::Point_2 */

/*!


Models the `ArrangementBasicTraits_2::X_monotone_curve_2` concept. 
Represents terminal segments of an algebraic curves, 
that means vertical segments or \f$ x\f$-monotone segments with no critical 
\f$ x\f$-coordinate in the interior of their \f$ x\f$-range. 
Terminal segments might either be bounded or unbounded. 
By definition, each interior point of 
a non-vertical segment has the same arc number (see the documentation of 
type `Point_2` above, which is called the <I>arc number</I> of the segment 
(note the arc number at the endpoints might differ). 
Such segments are represented internally by a 4-tuple \f$ (p,q,cv,arcno)\f$, 
where \f$ p\f$ and \f$ q\f$ are the endpoints, \f$ cv\f$ is the <I>supporting curve</I> 
that the segment belongs to, and arcno is the arc number of the segment. 

Arbitrary (weakly) \f$ x\f$-monotone segments are presented by a range 
of `X_monotone_curve_2` instances, whose union equals the segment. 
The functor `Construct_x_monotone_segment_2` allows their construction. 
To construct all (maximal) terminal segments of a curve, 
use the `Make_x_monotone_2` functor supplied by the traits class. 

*/
class Arr_algebraic_segment_traits_2::X_monotone_curve_2 {
public:

/// \name Modifiers 
/// @{

/*! 
returns the supporting algebraic curve of `s`. 
*/ 
Curve_2 curve () const; 

/*! 
returns whether `s` is a vertical segment. 
*/ 
bool is\_vertical () const; 

/*! 
returns whether `s` has a finite endpoint on the left 
*/ 
bool is\_finite (CGAL::Arr_curve_end ce) const; 

/*! 
\pre (The corresponding curve end is finite) 
*/ 
Point_2 curve_end (CGAL::Arr_curve_end ce) const; 

/*! 
returns the arc number of the segment. 
\pre (The segment is non-vertical) 
*/ 
int arcno () const; 

/*! 
returns the \f$ x\f$-coordinate of a vertical segment. 
\pre (The segment is vertical) 
*/ 
Algebraic_real_1 x () const; 

/// @}

}; /* end Arr_algebraic_segment_traits_2::X_monotone_curve_2 */

/*!


The `Curve_2` class nested within the B&eacute;zier traits class is used 
to represent a B&eacute;zier curve of arbitrary degree, which is defined by a 
sequence of rational control points. In addition to the methods listed 
below, the I/O operators `operator<<` and `operator>>` for 
standard output-streams are also supported. The copy constructor and 
assignment operator are supported as well. 

*/
class Arr_Bezier_curve_traits_2::Curve_2 {
public:

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
Curve_2 (); 

/*! 
constructs a B&eacute;zier curve as defined by the given range of control 
points. The value-type of `InputIterator` is `RatKernel::Point_2`. 
\pre The input range must contain at least two control points. 

*/ 
template <class InputIterator> 
Curve_2 (InputIterator pts_begin, InputIterator pts_end); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
returns the number of control points that define `B`. 
*/ 
size_t number_of_control_point () const; 

/*! 
returns the \f$ k\f$th control point. Note that the first control point equals 
the curve source, while the last control point equals its target. The rest 
of the control points do not lie on the curve. 
\pre \f$ k\f$ is smaller than the number of control points. 
*/ 
typename RatKernel::Point_2 control_point (size_t k) const; 

/*! 
returns the point \f$ B(t)\f$ on the curve that corresponds to the given 
rational parameter value. 
*/ 
typename RatKernel::Point_2 operator() (const Rational& t) const; 

/*! 
returns the point \f$ B(t)\f$ on the curve that corresponds to the given 
algebraic parameter value. 
*/ 
typename AlgKernel::Point_2 operator() (const Algebraic& t) const; 

/// @}

}; /* end Arr_Bezier_curve_traits_2::Curve_2 */

/*!


The `Point_2` class nested within the B&eacute;zier traits class is used 
to represent: (i) an endpoint of a B&eacute;zier curve, (ii) a vertical tangency 
point of a curve, used to subdivide it into \f$ x\f$-monotone subcurve, and 
(iii) an intersection point between two curves. While, points of type (i) have 
rational coordinates and are given as part of the input, points of the two 
latter types have algebraic coordinates. However, to speed up the arrangement 
construction, such point are not computed in an exact manner, and instead 
are given in an approximate representation. Note that the exact coordinates 
of a point may only be accessed if it is exactly computed. 

In addition to the methods listed below, the copy constructor and assignment 
operator for `Point_2` objects are also supported. 

*/
class Arr_Bezier_curve_traits_2::Point_2 {
public:

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
Point_2 (); 

/*! 
constructs the point \f$ B(t_0)\f$ on the given curve. As \f$ t_0\f$ is an 
algebraic number, the point has algebraic coordinates. 
*/ 
Point_2 (const Curve_2& B, const Algebraic& t_0); 

/*! 
constructs the point \f$ B(t_0)\f$ on the given curve. As \f$ t_0\f$ is a 
rational number, the point has rational coordinates. 
*/ 
Point_2 (const Curve_2& B, const Rational& t_0); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
returns the approximated coordinates of `p`. 
*/ 
std::pair<double, double> approximate () const; 

/*! 
returns whether the coordinates of `p` are computed in an exact manner. 
*/ 
bool is_exact () const; 

/*! 
returns the \f$ x\f$-coordinate of `p`. 
\pre `p` is exactly computed. 
*/ 
Algebraic x () const; 

/*! 
returns the \f$ y\f$-coordinate of `p`. 
\pre `p` is exactly computed. 
*/ 
Algebraic y () const; 

/*! 
returns whether the coordinates of `p` are rational numbers. 
*/ 
bool is_rational () const; 

/*! 
casts `p` to a point with rational coordinates. 
\pre `p` has rational coordinates. 
*/ 
operator typename RatKernel::Point_2 () const; 

/// @}

}; /* end Arr_Bezier_curve_traits_2::Point_2 */

/*!


The `X_monotone_curve_2` class nested within the B&eacute;zier traits is 
used to represent \f$ x\f$-monotone subcurves of B&eacute;zier curves. The subcurve is 
defined by a supporting B&eacute;zier curve \f$ B(t)\f$ and a range of definition in 
the parameter space \f$ [t_1, t_2] \subseteq [0, 1]\f$, where \f$ B(t_1)\f$ is the 
subcurve source and \f$ B(t_2)\f$ is its target. Note that as the point endpoints 
may only be approximated, the parameter range defining the subcurve may 
only be approximately known. 

It is not possible to construct \f$ x\f$-monotone subcurves directly. Instead, 
use the `Make_x_monotone_2` functor supplied by the traits class to 
subdivide a `Curve_2` object into \f$ x\f$-monotone subcurves. 

*/
class Arr_Bezier_curve_traits_2::X_monotone_curve_2 {
public:

/// \name Access Functions 
/// @{

/*! 
returns the supporting B&eacute;zier curve of `b`. 
*/ 
Curve_2 supporting_curve () const; 

/*! 
returns the source point of `b`. 
*/ 
Point_2 source () const; 

/*! 
returns the target point of `b`. 
*/ 
Point_2 target () const; 

/*! 
returns the left (\f$ xy\f$-lexicographically smaller) endpoint of `b`. 
*/ 
Point_2 left () const; 

/*! 
returns the right (\f$ xy\f$-lexicographically smaller) endpoint of `b`. 
*/ 
Point_2 right () const; 

/*! 
return the approximate parameter range defining the subcurve `b`. 
*/ 
std::pair<double, double> parameter_range () const; 

/// @}

}; /* end Arr_Bezier_curve_traits_2::X_monotone_curve_2 */

/*!


The `Curve_2` class nested within the traits class can represent 
arbitrary circular arcs, full circles and line segments and support their 
construction in various ways. 
The copy and default constructor as well as the assignment operator are 
provided. In addition, an `operator<<` for the curves is defined for 
standard output streams. 

*/
class Arr_circle_segment_traits_2::Curve_2 {
public:

/// \name Creation 
/// @{

/*! 
constructs an curve corresponding to the line segment `seg`. 
*/ 
Curve_2 (const typename Kernel::Segment_2& seg); 

/*! 
constructs an curve corresponding to the line segment directed 
from `source` to `target`, both having rational coordinates. 
*/ 
Curve_2 (const typename Kernel::Point_2& source, 
const typename Kernel::Point_2& target); 

/*! 
constructs an curve corresponding to the line segment supported by 
the given line, directed from `source` to `target`. 
Note that the two endpoints may have one-root coordinates. 
\pre Both endpoints must lie on the given supporting line. 
*/ 
Curve_2 (const typename Kernel::Line_2& line, 
const Point_2& source, 
const Point_2& target); 

/*! 
constructs an curve corresponding to the given circle. `circ` 
has a center point with rational coordinates and its <I>squared</I> 
radius is rational. 
*/ 
Curve_2 (const typename Kernel::Circle_2& circ); 

/*! 
constructs an curve corresponding to a circle centered at the rational 
point `c` whose radius `r` is rational. 
*/ 
Curve_2 (const typename Kernel::Point_2& c, 
const typename Kernel::FT& r, 
Orientation orient = COUNTERCLOCKWISE); 

/*! 
constructs a circular arc supported by `circ`, which has a 
center point with rational coordinates and whose <I>squared</I> 
radius is rational, with the given endpoints. The orientation of the 
arc is the same as the orientation of `circ`. 
\pre Both endpoints must lie on the given supporting circle. 
*/ 
Curve_2 (const typename Kernel::Circle_2& circ, 
const Point_2& source, const Point_2& target); 

/*! 
constructs a circular arc supported by a circle centered at the rational 
point `c` whose radius `r` is rational, directed from 
`source` to `target` with the given orientation. 
\pre Both endpoints must lie on the supporting circle. 
*/ 
Curve_2 (const typename Kernel::Point_2& c, 
const typename Kernel::FT& r, 
Orientation orient, 
const Point_2& source, const Point_2& target); 

/*! 
constructs an circular arc whose endpoints are `source` and 
`target` that passes through `mid`. All three points have 
rational coordinates. 
\pre The three points must not be collinear. 
*/ 
Curve_2 (const typename Kernel::Point_2& source, 
const typename Kernel::Point_2& mid, 
const typename Kernel::Point_2& target); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
indicates whether the curve represents a full circle. 
*/ 
bool is_full() const; 

/*! 
returns the source point. 
\pre `cv` is not a full circle. 
*/ 
const Point_2& source() const; 

/*! 
returns the target point. 
\pre `cv` is not a full circle. 
*/ 
const Point_2& target() const; 

/*! 
returns the orientation of the curve (`COLLINEAR` in case of 
line segments). 
*/ 
Orientation orientation() const; 

/*! 
determines whether `cv` is a line segment. 
*/ 
bool is_linear () const; 

/*! 
determines whether `cv` is a circular arc. 
*/ 
bool is_circular () const; 

/*! 
returns the supporting line of `cv`. 
\pre `cv` is a line segment. 
*/ 
typename Kernel::Line_2 supporting_line() const; 

/*! 
returns the supporting circle of `cv`. 
\pre `cv` is a circular arc. 
*/ 
typename Kernel::Circle_2 supporting_circle() const; 

/// @}

}; /* end Arr_circle_segment_traits_2::Curve_2 */

/*!


The `Point_2` number-type nested within the traits class represents 
a Cartesian point whose coordinates are algebraic numbers of type 
`CoordNT`. 

*/
class Arr_circle_segment_traits_2::Point_2 {
public:

/// \name Types 
/// @{

/*! 
the `Kernel::FT` type. 
*/ 
typedef Hidden_type Rational; 

/*! 
the algebraic number-type. 
*/ 
typedef Hidden_type CoordNT; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
Point_2 (); 

/*! 
creates the point \f$ (x,y)\f$. 
*/ 
Point_2 (const Rational& x, const Rational& y); 

/*! 
creates the point \f$ (x,y)\f$. 
*/ 
Point_2 (const CoordNT& x, const CoordNT& y); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
returns the \f$ x\f$-coordinate. 
*/ 
CoordNT x () const; 

/*! 
returns the \f$ y\f$-coordinate. 
*/ 
CoordNT y () const; 

/// @}

}; /* end Arr_circle_segment_traits_2::Point_2 */

/*!


The `X_monotone_curve_2` class nested within the traits class can 
represent \f$ x\f$-monotone and line segments (which are always weakly \f$ x\f$-monotone). 
The copy and default constructor as well as the assignment operator are 
provided. In addition, an `operator<<` for the curves is defined for 
standard output streams. 

*/
class Arr_circle_segment_traits_2::X_monotone_curve_2 {
public:

/// \name Creation 
/// @{

/*! 
constructs an curve corresponding to the line segment directed 
from `source` to `target`, both having rational coordinates. 
*/ 
X_monotone_curve_2 (const typename Kernel::Point_2& source, 
const typename Kernel::Point_2& target); 

/*! 
constructs an curve corresponding to the line segment supported by 
the given line, directed from `source` to `target`. 
Note that the two endpoints may have one-root coordinates. 
\pre Both endpoints must lie on the given supporting line. 
*/ 
X_monotone_curve_2 (const typename Kernel::Line_2& line, 
const Point_2& source, 
const Point_2& target); 

/*! 
constructs a circular arc supported by `circ`, which has a 
center point with rational coordinates and whose <I>squared</I> 
radius is rational, with the given endpoints. The orientation of the 
arc is determined by `orient`. 
\pre Both endpoints must lie on the given supporting circle. 
\pre The circular arc is \f$ x\f$-monotone. 
*/ 
X_monotone_curve_2 (const typename Kernel::Circle_2& circ, 
const Point_2& source, const Point_2& target, 
Orientation orient); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
returns the source point of `xcv`. 
*/ 
const Point_2& source() const; 

/*! 
returns the target point of `xcv`. 
*/ 
const Point_2& target() const; 

/*! 
returns true if `xcv` is directed right, false otherwise. 
*/ 
bool is_directed_right () const; 

/*! 
returns the left (lexicographically smaller) endpoint of `xcv`. 
*/ 
const Point_2& left() const; 

/*! 
returns the right (lexicographically larger) endpoint of `xcv`. 
*/ 
const Point_2& right() const; 

/*! 
returns the orientation of the curve (`COLLINEAR` in case of 
line segments). 
*/ 
Orientation orientation() const; 

/*! 
determines whether `xcv` is a line segment. 
*/ 
bool is_linear () const; 

/*! 
determines whether `xcv` is a circular arc. 
*/ 
bool is_circular () const; 

/*! 
returns the supporting line of `xcv`. 
\pre `xcv` is a line segment. 
*/ 
typename Kernel::Line_2 supporting_line() const; 

/*! 
returns the supporting circle of `xcv`. 
\pre `xcv` is a circular arc. 
*/ 
typename Kernel::Circle_2 supporting_circle() const; 

/*! 
returns a bounding box of the arc `xcv`. 
*/ 
Bbox_2 bbox() const; 

/// @}

}; /* end Arr_circle_segment_traits_2::X_monotone_curve_2 */

/*!


The `Curve_2` class nested within the conic-arc traits can represent 
arbitrary conic arcs and support their construction in various ways. 
The copy and default constructor as well as the assignment operator are 
provided for conic arcs. In addition, an `operator<<` 
for the curves is defined for standard output streams. 

*/
class Arr_conic_traits_2::Curve_2 {
public:

/// \name Creation 
/// @{

/*! 
constructs an arc corresponding to the line segment `seg`. 
*/ 
Curve_2 (const typename RatKernel::Segment_2& seg); 

/*! 
constructs an arc corresponding to the full circle `circ` 
(note that this circle has a center point with rational coordinates 
and rational squared radius). 
*/ 
Curve_2 (const typename RatKernel::Circle_2& circ); 

/*! 
constructs a circular arc supported by the circle `circ`, going 
in the given orientation `o` from the source point `ps` to 
its target point `pt`. 
\pre `ps` and `pt` both lie on the circle `circ`. 
\pre `o` is not `COLLINEAR`. 
*/ 
Curve_2 (const typename RatKernel::Circle_2& circ, 
Orientation o, 
const Point_2& ps, 
const Point_2& pt); 

/*! 
constructs a circular arc going from `p1` (its source point) 
through `p2` to `p3` (its target point). Note that all three 
points have rational coordinates. The orientation of the arc is 
determined automatically. 
\pre The three points are not collinear. 
*/ 
Curve_2 (const typename RatKernel::Point_2& p1, 
const typename RatKernel::Point_2& p2, 
const typename RatKernel::Point_2& p3); 

/*! 
constructs a conic arc that corresponds to the full conic curve 
\f$ r x^2 + s y^2 + t x y + u x + v y + w = 0\f$. 
\pre As a conic arc must be bounded, the given curve must be an ellipse, that is \f$ 4 r s - t^2 > 0\f$. 
*/ 
Curve_2 (const Rational& r, const Rational& s, 
const Rational& t, const Rational& u, 
const Rational& v, const Rational& w); 

/*! 
constructs a conic arc supported by the conic curve 
\f$ r x^2 + s y^2 + t x y + u x + v y + w = 0\f$, going 
in the given orientation `o` from the source point `ps` to 
its target point `pt`. 
\pre `ps` and `pt` both satisfy the equation of the supporting conic curve and define a bounded segment of this curve (e.g. in case of a hyperbolic arc, both point should be located on the same branch of the hyperbola). 
\pre `o` is not `COLLINEAR` if the supporting conic is curves, and must be `COLLINEAR` if it is not curved (a line or a line-pair). 
*/ 
Curve_2 (const Rational& r, const Rational& s, 
const Rational& t, const Rational& u, 
const Rational& v, const Rational& w, 
Orientation o, 
const Point_2& ps, 
const Point_2& pt); 

/*! 
constructs a conic arc going from `p1` (its source point) 
through `p2`, `p3` and `p4` (in this order) to `p5` 
(its target point). Note that all five points have rational coordinates. 
The orientation of the arc is determined automatically. 
\pre No three points of the five are not collinear. 
\pre The five points define a valid arc, in their given order. 
*/ 
Curve_2 (const typename RatKernel::Point_2& p1, 
const typename RatKernel::Point_2& p2, 
const typename RatKernel::Point_2& p3, 
const typename RatKernel::Point_2& p4, 
const typename RatKernel::Point_2& p5); 

/*! 
constructs a conic arc supported by the conic curve 
\f$ r x^2 + s y^2 + t x y + u x + v y + w = 0\f$, going 
in the given orientation `o` from its source point to its target 
Point. In this case only some approximations of the endpoints 
(`app_ps` and `app_pt`, respectively) is available, 
and their exact locations are given implicitly, specified by the 
intersections of the supporting conic curve with 
\f$ r_1 x^2 + s_1 y^2 + t_1 x y + u_1 x + v_1 y + w_1 = 0\f$ and 
\f$ r_2 x^2 + s_2 y^2 + t_2 x y + u_2 x + v_2 y + w_2 = 0\f$, respectively. 
\pre The two auxiliary curves specifying the endpoints really intersect with the supporting conic curve, such that the arc endpoints define a bounded segment of the supporting curve (e.g. in case of a hyperbolic arc, both point should be located on the same branch of the hyperbola). 
\pre `o` is not `COLLINEAR` if the supporting conic is curves, and must be `COLLINEAR` if it is not curved (a line or a line-pair). 
*/ 
Curve_2 (const Rational& r, const Rational& s, 
const Rational& t, const Rational& u, 
const Rational& v, const Rational& w, 
Orientation o, 
const Point_2& app_ps, 
const Rational& r1, const Rational& s1, 
const Rational& t1, const Rational& u1, 
const Rational& v1, const Rational& w1, 
const Point_2& app_pt, 
const Rational& r2, const Rational& s2, 
const Rational& t2, const Rational& u2, 
const Rational& v2, const Rational& w2); 

/// @} 

/// \name Access Functions 
CONVERROR Check if this needs to be spread\n/// The six following methods return the coefficients of the supported conic, after their conversion to integer number (represented by the `Integer` type of the `NtTraits` class):
/// @{

/*! 
indicates whether `a` is a valid conic arc. As the precondition to 
some of the constructor listed above are quite complicated, their 
violation does not cause the program to abort. Instead, the constructed 
arc is invalid (a defaultly constructed arc is also invalid). 
It is however recommended to check that a constructed arc is valid before 
inserting it to an arrangement, as this operation <I>will</I> cause the 
program to abort. 
*/ 
bool is_valid() const; 

/*! 
determines whether the arc is \f$ x\f$-monotone, namely each vertical line 
intersects it at most once. A vertical line segment is also considered 
(weakly) \f$ x\f$-monotone. 
*/ 
bool is_x_monotone() const; 

/*! 
determines whether the arc is \f$ y\f$-monotone, namely each horizontal line 
intersects it at most once. A horizontal line segment is also considered 
(weakly) \f$ x\f$-monotone. 
*/ 
bool is_y_monotone() const; 

/*! 
indicates whether the arc represents a full conic curve (en ellipse or 
a circle). 
*/ 
bool is_full_conic() const; 

/*! 
returns the coefficient of \f$ x^2\f$. 
*/ 
const typename NtTraits::Integer& r() const; 

/*! 
returns the coefficient of \f$ t^2\f$. 
*/ 
const typename NtTraits::Integer& s() const; 

/*! 
returns the coefficient of \f$ x y\f$. 
*/ 
const typename NtTraits::Integer& t() const; 

/*! 
returns the coefficient of \f$ x\f$. 
*/ 
const typename NtTraits::Integer& u() const; 

/*! 
returns the coefficient of \f$ y\f$. 
*/ 
const typename NtTraits::Integer& v() const; 

/*! 
returns the free coefficient. 
*/ 
const typename NtTraits::Integer& w() const; 

/*! 
returns the source point of the arc. 
\pre `a` is not a full conic curve. 
*/ 
const Point_2& source() const; 

/*! 
returns the target point of the arc. 
\pre `a` is not a full conic curve. 
*/ 
const Point_2& target() const; 

/*! 
returns the orientation of the arc. 
*/ 
Orientation orientation() const; 

/*! 
return a bounding box of the arc `a`. 
*/ 
Bbox_2 bbox() const; 

/// @} 

/// \name Operations 
/// @{

/*! 
sets a new source point for the conic arc. 
\pre `ps` lies on the supporting conic of `a`. 
*/ 
void set_source (const Point_2 & ps); 

/*! 
sets a new target point for the conic arc. 
\pre `pt` lies on the supporting conic of `a`. 
*/ 
void set_target (const Point_2 & pt); 

/// @}

}; /* end Arr_conic_traits_2::Curve_2 */

/*!


The `X_monotone_curve_2` class nested within the conic-arc traits is 
used to represent \f$ x\f$-monotone conic arcs. It inherits from the `Curve_2` 
type, therefore supports the access methods and the operations listed above. 

For efficiency reasons, we recommend users not to construct \f$ x\f$-monotone 
conic arc directly, but rather use the `Make_x_monotone_2` functor 
supplied by the conic-arc traits class to convert conic curves to 
\f$ x\f$-monotone curves. 

*/
class Arr_conic_traits_2::X_monotone_curve_2 {
public:

/// \name Creation 
/// @{

/*! 
converts the given arc to an \f$ x\f$-monotone arc. 
\pre `arc` is \f$ x\f$-monotone. 
*/ 
X_monotone_curve_2 (const Curve_2& arc); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
returns the left (lexicographically smaller) endpoint of `xa`. 
*/ 
const Point_2& left() const; 

/*! 
returns the right (lexicographically larger) endpoint of `xa`. 
*/ 
const Point_2& right() const; 

/// @}

}; /* end Arr_conic_traits_2::X_monotone_curve_2 */

/*!


The `Data_container` class nested within the consolidated 
curve-data traits and associated with the `Traits::X_monotone_curve_2` 
type is maintained as a list with unique data objects. This representation is 
simple and efficient in terms of memory consumption. It also requires that 
the `Data` class supports only the equality operator. Note however that 
most set operations require linear time. 

*/
class Arr_consolidated_curve_data_traits_2::Data_container {
public:

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
Data_container (); 

/*! 
constructs set containing a single `data` object. 
*/ 
Data_container (const Data& data); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
returns the number of data objects in the set. 
*/ 
std::size_t size () const; 

/*! 
returns an iterator pointing to the first data object. 
*/ 
Data_iterator begin () const; 

/*! 
returns a past-the-end iterator for the data objects. 
*/ 
Data_iterator end () const; 

/*! 
returns the first data object inserted into the set. 
\pre The number of data objects is not \f$ 0\f$. 
*/ 
const Data& front () const; 

/*! 
returns the last data object inserted into the set. 
\pre The number of data objects is not \f$ 0\f$. 
*/ 
const Data& back () const; 

/// @} 

/// \name Predicates 
/// @{

/*! 
check if the two sets contain the same data objects (regardless of order). 
*/ 
bool operator== (const Data_container& other) const; 

/*! 
find the given `data` object in the set and returns an iterator 
for this object, or `end()` if it is not found. 
*/ 
Data_iterator find (const Data& data); 

/// @} 

/// \name Modifiers 
/// @{

/*! 
inserts the given `data` object into the set. Returns `true` on 
success, or `false` if the set already contains the object. 
*/ 
bool insert (const Data& data); 

/*! 
erases the given `data` object from the set. Returns `true` on 
success, or `false` if the set does not contain the object. 
*/ 
bool erase (const Data& data); 

/*! 
clears the set. 
*/ 
void clear (); 

/// @}

}; /* end Arr_consolidated_curve_data_traits_2::Data_container */

/*!


The `Curve_2` class nested within the curve-data traits 
extends the `Base_traits_2::Curve_2` type with an extra data field of type 
`Data`. 

*/
class Arr_curve_data_traits_2::Curve_2
  : public Base_curve_2
 {
public:

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
Curve_2 (); 

/*! 
constructs curve from the given `base` curve with uninitialized 
data field. 
*/ 
Curve_2 (const Base_curve_2& base); 

/*! 
constructs curve from the given `base` curve with an attached 
`data` field. 
*/ 
Curve_2 (const Base_curve_2& base, const Data& data); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
returns the data field (a non-const version, which returns a reference 
to the data object, is also available). 
*/ 
const Curve_data& data () const; 

/*! 
sets the data field. 
*/ 
void set_data (const Curve_data& data); 

/// @}

}; /* end Arr_curve_data_traits_2::Curve_2 */

/*!


The `X_monotone_curve_2` class nested within the 
curve-data traits extends the `Base_traits_2::X_monotone_curve_2` type 
with an extra data field. 
*/
class Arr_curve_data_traits_2::X_monotone_curve_2 : public Base_x_monotone_curve_2 {
public:

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
X_monotone_curve_2 (); 

/*! 
constructs an \f$ x\f$-monotone curve from the given `base` curve with 
uninitialized data field. 
*/ 
X_monotone_curve_2 (const Base_x_monotone_curve_2& base); 

/*! 
constructs an \f$ x\f$-monotone curve from the given `base` \f$ x\f$-monotone 
curve with an attached `data` field. 
*/ 
X_monotone_curve_2 (const Base_x_monotone_curve_2& base, 
const X_monotone_curve_data& data); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
returns the field (a non-const version, which returns a reference 
to the data object, is also available). 
*/ 
const X_monotone_curve_data& data () const; 

/*! 
sets the data field. 
*/ 
void set_data (const X_monotone_curve_data& data); 

/// @}

}; /* end Arr_curve_data_traits_2::X_monotone_curve_2 */

/*!


The basic <span class="textsc">Dcel</span> face type. Serves as a basis class for an extended 
face record with auxiliary data fields. 

\models ::ArrangementDcelFace 

*/
class Arr_dcel_base::Arr_face_base {

}; /* end Arr_dcel_base::Arr_face_base */

/*!


The basic <span class="textsc">Dcel</span> halfedge type. Serves as a basis class for an 
extended halfedge record with auxiliary data fields. The `Curve` 
parameter is the type of \f$ x\f$-monotone curves associated with the vertices. 

\models ::ArrangementDcelHalfedge 

*/
template< typename Curve >
class Arr_dcel_base::Arr_halfedge_base {

}; /* end Arr_dcel_base::Arr_halfedge_base */

/*!


The basic <span class="textsc">Dcel</span> vertex type. Serves as a basis class for an extended 
vertex record with auxiliary data fields. The `Point` parameter is 
the type of points associated with the vertices. 

\models ::ArrangementDcelVertex 

*/
template< typename Point >
class Arr_dcel_base::Arr_vertex_base {

}; /* end Arr_dcel_base::Arr_vertex_base */

/*!


The `Curve_2` class nested within the polyline traits is used to 
represent general continuous piecewise-linear curves (a polyline can be 
self-intersecting) and support their construction from any range of points. 

The copy and default constructor as well as 
the assignment operator are provided for polyline curves. In addition, 
an `operator<<` for the curves is defined for standard output streams, 
and an `operator>>` for the curves is defined for standard input streams. 

*/
class Arr_polyline_traits_2::Curve_2 {
public:

/// \name Types 
/// @{

/*! 
A bidirectional iterator that allows 
traversing the points that comprise a polyline curve. 
*/ 
typedef Hidden_type const_iterator; 

/*! 
A bidirectional iterator that 
allows traversing the points that comprise a polyline curve. 
*/ 
typedef Hidden_type const_reverse_iterator; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor that constructs an empty polyline. 
*/ 
Curve_2 (); 

/*! 
constructs a polyline defined by the given range of points 
`[first, last)` (the value-type of `InputIterator` must be 
`SegmentTraits::Point_2`. 
If the range contains \f$ (n + 1)\f$ points labeled \f$ (p_{0},p_{1},\ldots,p_{n})\f$, 
the generated polyline consists of \f$ n\f$ segments, where the \f$ k\f$th segment 
is defined by the endpoints \f$ [p_{k-1},p_{k}]\f$. The first point in the 
range is considered as the source point of the polyline while the last 
point is considered as its target. 
\pre There are at least two points in the range. 
*/ 
template <class InputIterator> 
Curve_2 (Iterator first, Iterator last); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
returns the number of points that comprise the polyline. 
Note that if there are \f$ n\f$ points in the polyline, it is comprised 
of \f$ (n - 1)\f$ segments. 
*/ 
size_t points() const; 

/*! 
returns an iterator pointing at the source point of the polyline. 
*/ 
const_iterator begin() const; 

/*! 
returns an iterator pointing after the end of the polyline. 
*/ 
const_iterator end() const; 

/*! 
returns an iterator pointing at the target point of the polyline. 
*/ 
const_iterator rbegin() const; 

/*! 
returns an iterator pointing before the beginning of the polyline. 
*/ 
const_iterator rend() const; 

/*! 
returns the number of line segments comprising the polyline 
(equivalent to `pi.points() - 1`). 
*/ 
size_t size() const; 

/*! 
returns the \f$ k\f$th segment of the polyline. 
\pre `k` is not greater or equal to `pi.size() - 1`. 
*/ 
typename SegmentTraits::X_monotone_curve_2 
operator[] (size_t k) const; 

/*! 
return a bounding box of the polyline `pi`. 
*/ 
Bbox_2 bbox() const; 

/// @} 

/// \name Operations 
/// @{

/*! 
adds a new point to the polyline, which becomes the new target point 
of `pi`. 
*/ 
void push_back (const Point_2 & p); 

/*! 
resets the polyline. 
*/ 
void clear(); 

/// @}

}; /* end Arr_polyline_traits_2::Curve_2 */

/*!


The `X_monotone_curve_2` class nested within the polyline traits is used 
to represent \f$ x\f$-monotone piecewise linear curves. It inherits from the 
`Curve_2` type. It has a default constructor and a constructor from a 
range of points, just like the `Curve_2` class. However, there is 
precondition that the point range define an \f$ x\f$-monotone polyline. 

The points that define the \f$ x\f$-monotone polyline are 
always stored in an ascending lexicographical \f$ xy\f$-order, so their order may 
be reversed with respect to the input sequence. Also note that the 
\f$ x\f$-monotonicity ensures that an \f$ x\f$-monotone polyline is never 
self-intersecting (thus, a self-intersecting polyline will be subdivided 
to several interior-disjoint \f$ x\f$-monotone subcurves). 

*/
class Arr_polyline_traits_2::X_monotone_curve_2 {

}; /* end Arr_polyline_traits_2::X_monotone_curve_2 */

/*!


Functor to construct a `Curve_2`. To enable caching the class is not 
default constructible and must be obtained via the function 
`construct_curve_2_object()`, which is a member of the traits. 

\models ::Assignable 
\models ::CopyConstructible 
\models ::AdaptableBinaryFunction 
\models ::AdaptableUnaryFunction 

*/
class Arr_rational_function_traits_2::Construct_curve_2 {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 Polynomial_1; 

/*! 

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 Algebraic_real_1; 

/*! 

*/ 
typedef Arr_rational_function_traits_2<AlgebraicKernel_d_1>::Curve_2 result_type; 

/*! 

*/ 
typedef Polynomial_1 argument_type; 

/*! 

*/ 
typedef Polynomial_1 first_argument_type; 

/*! 

*/ 
typedef Polynomial_1 second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*! 
Constructs a curve representing the polynomial function \f$ y = P(x)\f$. 
*/ 
Curve_2 operator()(Polynomial_1 P) const; 

/*! 
Constructs a curve representing the polynomial function \f$ y = P(x)\f$. 
The function is defined over the interval \f$ [x,+\infty)\f$ if \f$ right\f$ is true 
and \f$ (-\infty,x]\f$ otherwise. 
*/ 
Curve_2 operator()(Polynomial_1 P, const Algebraic_real_1& x, bool right) const; 

/*! 
Constructs a curve representing the polynomial function \f$ y = P(x)\f$. 
The function is defined over the interval \f$ [lower,upper]\f$. 
*/ 
Curve_2 operator()(Polynomial_1 P, const Algebraic_real_1& lower, const Algebraic_real_1& upper) const; 

/*! 
Constructs a curve representing the rational function \f$ y = P(x)/Q(x)\f$. 
*/ 
Curve_2 operator()(Polynomial_1 P, Polynomial_1 Q) const; 

/*! 
Constructs a curve representing the rational function \f$ y = P(x)/Q(x)\f$. 
The function is defined over the interval \f$ I=[x,+\infty)\f$ if \f$ right\f$ is 
true and \f$ I=(-\infty,x]\f$ otherwise. 
*/ 
Curve_2 operator()(Polynomial_1 P, Polynomial_1 Q, const Algebraic_real_1& x, bool right) const; 

/*! 
Constructs a curve representing the rational function \f$ y = P(x)/Q(x)\f$. 
The function is defined over the interval \f$ I=[lower,upper]\f$. 
*/ 
Curve_2 operator()(Polynomial_1 P, Polynomial_1 Q, const Algebraic_real_1& lower, const Algebraic_real_1& upper) const; 

/*! 
Constructs a curve representing the polynomial function \f$ y = P(x)\f$, where 
the coefficients of \f$ P\f$ are given in the range `[begin,end)`. 
*/ 
template <typename InputIterator> 
Curve_2 operator()(InputIterator begin, InputIterator end) const; 

/*! 
Constructs a curve representing the polynomial function \f$ y = P(x)\f$, where 
the coefficients of \f$ P\f$ are given in the range `[begin,end)`. The 
function is defined over the interval \f$ [x,+\infty)\f$ if \f$ right\f$ is true 
and \f$ (-\infty,x]\f$ otherwise. 
*/ 
template <typename InputIterator> 
Curve_2 operator()(InputIterator begin, InputIterator end, 
const Algebraic_real_1& x, bool right) const; 

/*! 
Constructs a curve representing the polynomial function \f$ y = P(x)\f$, where 
the coefficients of \f$ P\f$ are given in the range `[begin,end)`. The 
function is defined over the interval \f$ [lower,upper]\f$. 
*/ 
template <typename InputIterator> 
Curve_2 operator()(InputIterator begin, InputIterator end, 
const Algebraic_real_1& lower, 
const Algebraic_real_1& upper) const; 

/*! 
Constructs a curve representing the rational function \f$ y = P(x)/Q(x)\f$, 
where the coefficients of \f$ P\f$ and \f$ Q\f$ are given in the ranges 
`[begin_numer,end_numer)` and `[begin_denom,end_denom)`, 
respectively. 
*/ 
template <typename InputIterator> 
Curve_2 operator()(InputIterator begin_numer, InputIterator end_numer, 
InputIterator begin_denom, InputIterator end_denom) const; 

/*! 
Constructs a curve representing the rational function \f$ y = P(x)/Q(x)\f$, 
where the coefficients of \f$ P\f$ and \f$ Q\f$ are given in the ranges 
`[begin_numer,end_numer)` and `[begin_denom,end_denom)`, 
respectively. The function is defined over the interval \f$ I=[x,+\infty)\f$ 
if \f$ right\f$ is true and \f$ I=(-\infty,x]\f$ otherwise. 
*/ 
template <typename InputIterator> 
Curve_2 operator()(InputIterator begin_numer, InputIterator end_numer, 
InputIterator begin_denom, InputIterator end_denom, 
const Algebraic_real_1& x, bool right) const; 

/*! 
Constructs a curve representing the rational function \f$ y = P(x)/Q(x)\f$, 
where the coefficients of \f$ P\f$ and \f$ Q\f$ are given in the ranges 
`[begin_numer,end_numer)` and `[begin_denom,end_denom)`, 
respectively. The function is defined over the interval \f$ I=[lower,upper]\f$. 
*/ 
template <typename InputIterator> 
Curve_2 operator()(InputIterator begin_numer, InputIterator end_numer, 
InputIterator begin_denom, InputIterator end_denom, 
const Algebraic_real_1& lower, 
const Algebraic_real_1& upper) const; 

/// @}

}; /* end Arr_rational_function_traits_2::Construct_curve_2 */

/*!


Functor to construct a `X_monotone_curve_2`. To enable caching the class 
is not default constructible and must be obtained via the function 
`construct_x_monotone_curve_2_object()`, which is a member of the traits. 

\models ::Assignable 
\models ::CopyConstructible 
\models ::AdaptableBinaryFunction 
\models ::AdaptableUnaryFunction 

*/
class Arr_rational_function_traits_2::Construct_x_monotone_curve_2 {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 Polynomial_1; 

/*! 

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 Algebraic_real_1; 

/*! 

*/ 
typedef Arr_rational_function_traits_2<AlgebraicKernel_d_1>::X_monotone_curve_2 result_type; 

/*! 

*/ 
typedef Polynomial_1 argument_type; 

/*! 

*/ 
typedef Polynomial_1 first_argument_type; 

/*! 

*/ 
typedef Polynomial_1 second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*! 
Constructs an \f$ x\f$-monotone curve supported by the polynomial function 
\f$ y = P(x)\f$. 
*/ 
X_monotone_curve_2 operator()(Polynomial_1 P) const; 

/*! 
Constructs an \f$ x\f$-monotone curve supported by the polynomial function 
\f$ y = P(x)\f$. The function is defined over the interval \f$ [x,+\infty)\f$ if 
\f$ right\f$ is true and \f$ (-\infty,x]\f$ otherwise. 
*/ 
X_monotone_curve_2 operator()(Polynomial_1 P, 
const Algebraic_real_1& x, 
bool right) const; 

/*! 
Constructs an \f$ x\f$-monotone curve supported by the polynomial function 
\f$ y = P(x)\f$. The function is defined over the interval \f$ [lower,upper]\f$. 
*/ 
X_monotone_curve_2 operator()(Polynomial_1 P, 
const Algebraic_real_1& lower, 
const Algebraic_real_1& upper); const 

/*! 
Constructs an \f$ x\f$-monotone curve supported by the rational function 
\f$ y = P(x)/Q(x)\f$. 
\pre \f$ Q\f$ has no real roots. 
*/ 
X_monotone_curve_2 operator()(Polynomial_1 P, Polynomial_1 Q); const 

/*! 
Constructs an \f$ x\f$-monotone curve supported by the rational function 
\f$ y = P(x)/Q(x)\f$. The function is defined over the interval \f$ I=[x,+\infty)\f$ 
if \f$ right\f$ is true and \f$ I=(-\infty,x]\f$ otherwise. 
\pre \f$ Q\f$ has no real roots in the interior of \f$ I\f$. 
*/ 
X_monotone_curve_2 operator()(Polynomial_1 P, Polynomial_1 Q, 
const Algebraic_real_1& x, 
bool right); const 

/*! 
Constructs an \f$ x\f$-monotone curve supported by the rational function 
\f$ y = P(x)/Q(x)\f$. The function is defined over the interval \f$ I=[lower,upper]\f$. 
\pre \f$ Q\f$ has no real roots in the interior of \f$ I\f$. 
*/ 
X_monotone_curve_2 operator()(Polynomial_1 P, Polynomial_1 Q, 
const Algebraic_real_1& lower, 
const Algebraic_real_1& upper); const 

/*! 
Constructs an \f$ x\f$-monotone curve supported by the polynomial function 
\f$ y = P(x)\f$, where the coefficients of \f$ P\f$ are given in the range 
`[begin,end)`. 
*/ 
template <typename InputIterator> 
X_monotone_curve_2 operator()(InputIterator begin, InputIterator end) const; 

/*! 
Constructs an \f$ x\f$-monotone curve supported by the polynomial function 
\f$ y = P(x)\f$, where the coefficients of \f$ P\f$ are given in the range 
`[begin,end)`. The function is defined over the interval \f$ [x,+\infty)\f$ 
if \f$ right\f$ is true and \f$ (-\infty,x]\f$ otherwise. 
*/ 
template <typename InputIterator> 
X_monotone_curve_2 operator()(InputIterator begin, InputIterator end, 
const Algebraic_real_1& x, bool right) const; 

/*! 
Constructs an \f$ x\f$-monotone curve supported by the polynomial function 
\f$ y = P(x)\f$, where the coefficients of \f$ P\f$ are given in the range 
`[begin,end)`. The function is defined over the interval 
\f$ [lower,upper]\f$. 
*/ 
template <typename InputIterator> 
X_monotone_curve_2 operator()(InputIterator begin, InputIterator end 
const Algebraic_real_1& lower, 
const Algebraic_real_1& upper); const 

/*! 
Constructs an \f$ x\f$-monotone curve supported by the rational function 
\f$ y = P(x)/Q(x)\f$, where the coefficients of \f$ P\f$ and \f$ Q\f$ are given in the 
ranges `[begin_numer,end_numer)` and `[begin_denom,end_denom)`, 
respectively. 
\pre \f$ Q\f$ has no real roots. 
*/ 
template <typename InputIterator> 
X_monotone_curve_2 operator()(InputIterator begin_numer, InputIterator end_numer, 
InputIterator begin_denom, InputIterator end_denom); const 

/*! 
Constructs an \f$ x\f$-monotone curve supported by the rational function 
\f$ y = P(x)/Q(x)\f$, where the coefficients of \f$ P\f$ and \f$ Q\f$ are given in the 
ranges `[begin_numer,end_numer)` and `[begin_denom,end_denom)`, 
respectively. The function is defined over the interval \f$ I=[x,+\infty)\f$ 
if \f$ right\f$ is true and \f$ I=(-\infty,x]\f$ otherwise. 
\pre \f$ Q\f$ has no real roots in the interior of \f$ I\f$. 
*/ 
template <typename InputIterator> 
X_monotone_curve_2 operator()(InputIterator begin_numer, InputIterator end_numer, 
InputIterator begin_denom, InputIterator end_denom, 
const Algebraic_real_1& x, bool right); const 

/*! 
Constructs an \f$ x\f$-monotone curve supported by the rational function 
\f$ y = P(x)/Q(x)\f$, where the coefficients of \f$ P\f$ and \f$ Q\f$ are given in the 
ranges `[begin_numer,end_numer)` and `[begin_denom,end_denom)`, 
respectively. The function is defined over the interval \f$ I=[lower,upper]\f$. 
\pre \f$ Q\f$ has no real roots in the interior of \f$ I\f$. 
*/ 
template <typename InputIterator> 
X_monotone_curve_2 operator()(InputIterator begin_numer, InputIterator end_numer, 
InputIterator begin_denom, InputIterator end_denom, 
const Algebraic_real_1& lower, const Algebraic_real_1& upper); const 

/// @}

}; /* end Arr_rational_function_traits_2::Construct_x_monotone_curve_2 */

/*!


The `Curve_2` class nested within the traits is used 
to represent rational functions which may be restricted to a certain x-range. 

\models ::ArrTraits::Curve_2 

*/
class Arr_rational_function_traits_2::Curve_2 {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 Polynomial_1; 

/*! 

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 Algebraic_real_1; 

/// @} 

/// \name Operations 
/// @{

/*! 
returns the numerator of the supporting rational function. 
*/ 
const Polynomial_1& numerator () const; 

/*! 
returns the denominator of the supporting rational function. 
*/ 
const Polynomial_1& denominator () const; 

/*! 
returns whether &ccedil;urve is continuous, namely whether it does not 
contains any vertical asymptotes in its interior. 
*/ 
bool is_continuous() const; 

/*! 
returns whether the \f$ x\f$-coordinate of &ccedil;urve's left end is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space left_parameter_space_in_x () const; 

/*! 
returns whether the \f$ x\f$-coordinate of &ccedil;urve's right end is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space right_parameter_space_in_x () const; 

/*! 
returns the \f$ x\f$-coordinate of the left end. 
\pre left_boundary_in_x()==ARR_INTERIOR 
*/ 
Algebraic_real_1 left_x() const; 

/*! 
returns the \f$ x\f$-coordinate of the right end. 
\pre right_boundary_in_x()==ARR_INTERIOR 
*/ 
Algebraic_real_1 right_x() const; 

/// @}

}; /* end Arr_rational_function_traits_2::Curve_2 */

/*!


\models ::ArrTraits::Point_2 

*/
class Arr_rational_function_traits_2::Point_2 {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 Polynomial_1; 

/*! 

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 Algebraic_real_1; 

/*! 

*/ 
typedef AlgebraicKernel_d_1::Bound Bound; 

/// @} 

/// \name Operations 
/// @{

/*! 
returns the numerator of the supporting rational function. 
*/ 
Polynomial_1 numerator () const; 

/*! 
returns the denominator of the supporting rational function. 
*/ 
Polynomial_1 denominator () const; 

/*! 
returns double-approximations of the x- and y-coordinates. 
*/ 
std::pair<double,double> to_double() const; 

/*! 
returns the \f$ x\f$-coordinate of the point. 
*/ 
Algebraic_real_1 x() const; 

/*! 
obtains the y-coordinates of the point. <B>Attention:</B> As described above, 
points are not stored by their y-coordinate in `Algebraic_real_1` 
representation. In fact, this representation must be computed on demand, and 
might become quite costly for points defined by high-degree polynomials. 
Therefore, it is recommended to avoid calls to this function as much as 
possible. 
*/ 
Algebraic_real_1 y() const; 

/*! 
Computes a pair \f$ p\f$ approximating the \f$ x\f$-coordinate with 
respect to the given absolute precision \f$ a\f$. 
\post \f$ p.first \leq x \leq p.second \f$ 
\post \f$ p.second - p.first \leq2^{-a} \f$ 
*/ 
std::pair<Bound,Bound> approximate_absolute_x(int a) const; 

/*! 
Computes a pair \f$ p\f$ approximating the \f$ y\f$-coordinate with 
respect to the given absolute precision \f$ a\f$. 
\post \f$ p.first \leq y \leq p.second \f$ 
\post \f$ p.second - p.first \leq2^{-a} \f$ 
*/ 
std::pair<Bound,Bound> approximate_absolute_y(int a) const; 

/*! 
Computes a pair \f$ p\f$ approximating the \f$ x\f$-coordinate with 
respect to the given relative precision \f$ r\f$. 
\post \f$ p.first \leq x \leq p.second \f$ 
\post \f$ p.second - p.first \leq2^{-r}|x| \f$ 
*/ 
std::pair<Bound,Bound> approximate_relative_x(int r) const; 

/*! 
Computes a pair \f$ p\f$ approximating the \f$ y\f$-coordinate with 
respect to the given relative precision \f$ r\f$. 
\post \f$ p.first \leq y \leq p.second \f$ 
\post \f$ p.second - p.first \leq2^{-r}|y| \f$ 
*/ 
std::pair<Bound,Bound> approximate_relative_y(int r) const; 

/// @}

}; /* end Arr_rational_function_traits_2::Point_2 */

/*!


The `X_monotone_curve_2` class nested within the traits is used 
to represent \f$ x\f$-monotone parts of rational functions. In particular, such an \f$ x\f$-monotone curve 
may not contain a vertical asymptote in its interior \f$ x\f$-range. 

\models ::ArrTraits::XMonotoneCurve_2 

*/
class Arr_rational_function_traits_2::X_monotone_curve_2 {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef AlgebraicKernel_d_1::Polynomial_1 Polynomial_1; 

/*! 

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 Algebraic_real_1; 

/*! 

*/ 
typedef Arr_rational_function_traits_2<AlgebraicKernel_d_1>::Point_2 Point_2; 

/// @} 

/// \name Operations 
/// @{

/*! 
returns the numerator of the supporting rational function. 
*/ 
const Polynomial_1& numerator () const; 

/*! 
returns the denominator of the supporting rational function. 
*/ 
const Polynomial_1& denominator () const; 

/*! 
returns whether the \f$ x\f$-coordinate of the source is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space source_parameter_space_in_x () const; 

/*! 
returns whether the \f$ y\f$-coordinate of the source is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space source_parameter_space_in_y () const; 

/*! 
returns the source point of the arc. 
\pre Both the \f$ x\f$- and \f$ y\f$-coordinates of the source point is finite. 
*/ 
const Point_2& source() const; 

/*! 
returns the \f$ x\f$-coordinate of the source point. 
\pre The \f$ x\f$-coordinate of the source point is finite. 
*/ 
Algebraic_real_1 source_x() const; 

/*! 
returns whether the \f$ x\f$-coordinate of the target is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space target_parameter_space_in_x () const; 

/*! 
returns whether the \f$ y\f$-coordinate of the target is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space target_parameter_space_in_y () const; 

/*! 
returns the target point of the arc. 
\pre Both the \f$ x\f$- and \f$ y\f$-coordinates of the target point is finite. 
*/ 
const Point_2& target() const; 

/*! 
returns the \f$ x\f$-coordinate of the target point. 
\pre The \f$ x\f$-coordinate of the target point is finite. 
*/ 
Algebraic_real_1 target_x() const; 

/*! 
returns whether the \f$ x\f$-coordinate of the left curve end is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space left_parameter_space_in_x () const; 

/*! 
returns whether the \f$ y\f$-coordinate of the left curve end is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space left_parameter_space_in_y () const; 

/*! 
returns the left point of the arc. 
\pre Both the \f$ x\f$- and \f$ y\f$-coordinates of the left point is finite. 
*/ 
const Point_2& left() const; 

/*! 
returns the \f$ x\f$-coordinate of the left point. 
\pre The \f$ x\f$-coordinate of the left point is finite. 
*/ 
Algebraic_real_1 left_x() const; 

/*! 
returns whether the \f$ x\f$-coordinate of the right curve end is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space right_parameter_space_in_x () const; 

/*! 
returns whether the \f$ y\f$-coordinate of the right curve end is finite or 
whether it is \f$ \pm\infty\f$. 
*/ 
Arr_parameter_space right_parameter_space_in_y () const; 

/*! 
returns the right point of the arc. 
\pre Both the \f$ x\f$- and \f$ y\f$-coordinates of The right point is finite. 
*/ 
const Point_2& right() const; 

/*! 
returns the \f$ x\f$-coordinate of the right point. 
\pre The \f$ x\f$-coordinate of the right point is finite. 
*/ 
Algebraic_real_1 right_x() const; 

/*! 
returns whether the curve is oriented from left to right. 
*/ 
bool is_left_to_right () const; 

/// @}

}; /* end Arr_rational_function_traits_2::X_monotone_curve_2 */
