
namespace CGAL {

/*!
\ingroup PkgArrangement2

The class `Arr_conic_traits_2` is a model of the `ArrangementTraits_2` concept 
and can be used to construct and maintain arrangements of bounded segments of 
algebraic curves of degree \f$ 2\f$ at most, also known as <I>conic curves</I>. 

A general conic curve \f$ C\f$ is the locus of all points \f$ (x,y)\f$ satisfying the 
equation: \f$ r x^2 + s y^2 + t x y + u x + v y + w = 0\f$, where: 
<UL> 
<LI>If \f$ 4 r s - t^2 > 0\f$, \f$ C\f$ is an ellipse. 
A special case occurs when \f$ r = s\f$ and \f$ t = 0\f$, when \f$ C\f$ 
becomes a circle. 
<LI>If \f$ 4 r s - t^2 < 0\f$, \f$ C\f$ is a hyperbola. 
<LI>If \f$ 4 r s - t^2 = 0\f$, \f$ C\f$ is a parabola. 
A degenerate case occurs when \f$ r = s = t = 0\f$, when \f$ C\f$ is a line. 
</UL> 

A <I>bounded conic arc</I> is defined as either of the following: 
<UL> 
<LI>A full ellipse (or a circle) \f$ C\f$. 
<LI>The tuple \f$ \langle C, p_s, p_t, o \rangle\f$, where \f$ C\f$ is the supporting 
conic curve, with the arc endpoints being \f$ p_s\f$ and \f$ p_t\f$ 
(the source and target points, respectively). The orientation \f$ o\f$ 
indicates whether we proceed from \f$ p_s\f$ to \f$ p_t\f$ in a clockwise or in a 
counterclockwise direction. Note that \f$ C\f$ may also 
correspond to a line or to pair of lines - in this case \f$ o\f$ may 
specify a `COLLINEAR` orientation. 
</UL> 

A very useful subset of the set of conic arcs are line segments and circular 
arcs, as arrangements of circular arcs and line segments have some 
interesting applications (e.g. offsetting polygons, motion planning for a 
disc robot, etc.). Circular arcs and line segment are simpler objects and can 
be dealt with more efficiently than arbitrary arcs. For these reasons, it is 
possible to construct conic arcs from segments and from circles. Using these 
constructors is highly recommended: It is more straightforward and also speeds 
up the arrangement construction. However, in case the set of input curves 
contain only circular arcs and line segments, it is recommended to use the 
`Arr_circle_segment_2` class to achieve faster running times. 

In our representation, all conic coefficients (namely \f$ r, s, t, u, v, w\f$) 
must be rational numbers. This guarantees that the coordinates of all 
arrangement vertices (in particular, those representing intersection 
points) are algebraic numbers of degree \f$ 4\f$ (a real number \f$ \alpha\f$ 
is an algebraic number of degree \f$ d\f$ if there exist a polynomial \f$ p\f$ with 
<I>integer</I> coefficient of degree \f$ d\f$ such that \f$ p(\alpha) = 0\f$). 
We therefore require separate representations of the curve coefficients and 
the point coordinates. The `NtTraits` should be instantiated with a class 
that defines nested `Integer`, `Rational` and `Algebraic` 
number types and supports various operations on them, yielding certified 
computation results (for example, it can convert rational numbers to algebraic 
numbers and can compute roots of polynomials with integer coefficients). 
The other template parameters, `RatKernel` and `AlgKernel` should be 
geometric kernels templated with the `NtTraits::Rational` and 
`NtTraits::Algebraic` number types, respectively. 
It is recommended to instantiate the `CORE_algebraic_number_traits` 
class as the `NtTraits` parameter, with 
`Cartesian<NtTraits::Rational>` and `Cartesian<NtTraits::Algebraic>` 
instantiating the two kernel types, respectively. 
The number types in this case are provided by the <span class="textsc">Core</span> library, with its 
ability to exactly represent simple algebraic numbers. 

The traits class inherits its point type from `AlgKernel::Point_2`, 
and defines a curve and \f$ x\f$-monotone curve types, as detailed below. 

While the `Arr_conic_traits_2` models the concept 
`ArrangementDirectionalXMonotoneTraits_2`, the implementation of 
the `Arr_mergeable_2` operation does not enforce the input curves 
to have the same direction as a precondition. Moreover, `Arr_conic_traits_2` 
supports the merging of curves of opposite directions. 

\models ::ArrangementTraits_2 
\models ::ArrangementLandmarkTraits_2 
\models ::ArrangementDirectionalXMonotoneTraits_2 

Types 
-------------- 

CONVERROR 2 nested classes missing 

*/
template< typename RatKernel, typename AlgKernel, typename NtTraits >
class Arr_conic_traits_2 {
public:

/// \name Types 
/// @{

/*! 
the `NtTraits::Rational` type 
(and also the `RatKernel::FT` type). 
*/ 
typedef Hidden_type Rational; 

/*! 
the `NtTraits::Algebraic` type 
(and also the `AlgKernel::FT` type). 
*/ 
typedef Hidden_type Algebraic; 

/// @}

}; /* end Arr_conic_traits_2 */
} /* end namespace CGAL */
