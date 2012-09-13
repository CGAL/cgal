
namespace CGAL {

/*!
\ingroup PkgArrangement2TraitsClasses

The traits class `Arr_polyline_traits_2` is a model of the `ArrangementTraits_2` 
concept. It handles piecewise linear curves, commonly referred to as 
polylines. Each polyline is a chain of segments, where each two neighboring 
segments in the chain share a common endpoint. The traits class exploits the 
functionality of the `SegmentTraits` template-parameter to handle the 
segments that comprise the polyline curves. 

The class instantiated for the template parameter `SegmentTraits` must 
be a model of the `ArrangementTraits_2` concept that handles line 
segments (e.g., `Arr_segment_traits_2<Kernel>` or 
`Arr_non_caching_segment_traits_2<Kernel>`, where the first 
alternative is recommended). 

The number type used by the injected segment traits should support exact 
rational arithmetic (that is, the number type should support 
the arithmetic operations \f$ +\f$, \f$ -\f$, \f$ \times\f$ and \f$ \div\f$ that should be 
carried out without loss of precision), in order to avoid robustness 
problems, although other inexact number types could be used at the user's 
own risk. 

\models ::ArrangementTraits_2 
\models ::ArrangementLandmarkTraits_2 

\sa `Arr_segment_traits_2<Kernel>` 
\sa `Arr_non_caching_segment_traits_2<Kernel>` 

*/
template< typename SegmentTraits >
class Arr_polyline_traits_2 {
public:


/*!
The `Curve_2` class nested within the polyline traits is used to 
represent general continuous piecewise-linear curves (a polyline can be 
self-intersecting) and support their construction from any range of points. 

The copy and default constructor as well as 
the assignment operator are provided for polyline curves. In addition, 
an `operator<<` for the curves is defined for standard output streams, 
and an `operator>>` for the curves is defined for standard input streams. 

*/
class Curve_2 {
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
class X_monotone_curve_2 {

}; /* end Arr_polyline_traits_2::X_monotone_curve_2 */

}; /* end Arr_polyline_traits_2 */
} /* end namespace CGAL */
