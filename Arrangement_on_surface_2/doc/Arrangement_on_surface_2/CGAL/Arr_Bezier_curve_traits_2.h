
namespace CGAL {

/*!
\ingroup PkgArrangement2

The traits class `Arr_Bezier_curve_traits_2` is a model of the `ArrangementTraits_2` 
concept that handles planar B&eacute;zier curves. A planar <I>B&eacute;zier curve</I> 
\f$ B\f$ is a parametric curve defined by a sequence of <I>control points</I> 
\f$ p_0, \ldots, p_n\f$ as follows: 


\f{eqnarray*}{
B(t) = \left(X(t), Y(t)\right) 
= \ccSum{k=0}{n}{p_k \cdot \frac{n!}{k! (n-k)!} \cdot 
t^k (1-t)^{n-k}}\ . 
\f}
where \f$ t \in [0, 1]\f$. The degree of the curve is therefore \f$ n\f$ - 
namely, \f$ X(t)\f$ and \f$ Y(t)\f$ are polynomials of degree \f$ n\f$. B&eacute;zier curves 
have numerous applications in computer graphics and solid modelling. They 
are used, for example, in free-form sketches and for defining the true-type 
fonts. 

In our representation, we assume that the coordinates of all control 
points are rational numbers (namely they are given as objects of the 
`RatKernel::Point_2` type), so both \f$ X(t)\f$ and \f$ Y(t)\f$ are polynomials 
with rational coefficients. The intersection points between curves are 
however algebraic numbers, and their exact computation is time-consuming. 
The traits class therefore contains a layer of geometric filtering that 
performs all computation in an approximate manner whenever possible, and 
it resorts to exact computations only when the approximate computation 
fails to produce an unambiguous result. 

We therefore require separate representations of the control points and 
the intersection points. The `NtTraits` should be instantiated with a class 
that defines nested `Integer`, `Rational` and `Algebraic` number 
types and supports various operations on them, yielding certified computation 
results (for example, in can convert rational numbers to algebraic numbers 
and can compute roots of polynomials with integer coefficients). 
The other template parameters, `RatKernel` and `AlgKernel` should be 
geometric kernels templated with the `NtTraits::Rational` and 
`NtTraits::Algebraic` number types, respectively. It is recommended to 
instantiate the `CORE_algebraic_number_traits` class as the `NtTraits` 
parameter, with `Cartesian<NtTraits::Rational>` and 
`Cartesian<NtTraits::Algebraic>` instantiating the two kernel types, 
respectively. The number types in this case are provided by the <span class="textsc">Core</span> 
library, with its ability to exactly represent simple algebraic numbers. 

While `Arr_Bezier_curve_traits_2` models the concept 
`ArrangementDirectionalXMonotoneTraits_2`, the implementation of 
the `Arr_mergeable_2` operation does not enforce the input curves 
to have the same direction as a precondition. Moreover, `Arr_Bezier_curve_traits_2` 
supports the merging of curves of opposite directions. 

\models ::ArrangementTraits_2 
\models ::ArrangementDirectionalXMonotoneTraits_2 

CONVERROR 3 nested classes missing 

*/
template< typename RatKernel, typename AlgKernel, typename NtTraits >
class Arr_Bezier_curve_traits_2 {
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

}; /* end Arr_Bezier_curve_traits_2 */
} /* end namespace CGAL */
