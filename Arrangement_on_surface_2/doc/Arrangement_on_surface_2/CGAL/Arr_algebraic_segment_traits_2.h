
namespace CGAL {

/*!
\ingroup PkgArrangement2

The traits class `Arr_algebraic_segment_traits_2` is a model of the `ArrangementTraits_2` 
concept that handles planar algebraic curves of arbitrary degree, 
and \f$ x\f$-monotone of such curves. 
A planar (real) <I>algebraic curve</I> 
is the vanishing set of a polynomial in two variables, that is, the 
curve is defined by the defining equation 
\f[ f(x):=\sum_{i+j\leq n} a_{ij} x^i y^j =0, \f] 
where \f$ n\f$ is the degree of the curve. 

The traits class allows the construction of algebraic curves, 
by specifying their implicit equation. \f$ x\f$-monotone and vertical segments 
of a curve can also be defined; unbounded curves and segments are supported. 
The template parameter `Coefficient` defines 
the innermost coefficient type of the polynomials. Currently, 
the types `leda::integer` and `CORE::BigInt` are supported as well 
as any instance of `CGAL::Sqrt_extension` that is instantiated with 
one of the integral types above. 

\models ::ArrangementTraits_2 

6 nested classes missing

*/
template< typename Coefficient >
class Arr_algebraic_segment_traits_2 {
public:

/// \name Types 
/// @{

/*! 
Value to specify whether a point should be in the interior 
of a segment, or its minimal point, 
or its maximal point in lexicographic order. 
*/
enum Site_of_point { POINT_IN_INTERIOR = 0, MIN_ENDPOINT = -1, MAX_ENDPOINT = 1 };

/*! 
the type for bivariate polynomials, 
with innermost coefficient type `Coefficient`. 
Constitutes a model of the concept `Polynomial_d` 
with two variables. 

\sa CGAL::Polynomial_d
*/ 
typedef Hidden_type Polynomial_2; 

/*! 
model for the concept 
`AlgebraicKernel_1` 
*/ 
typedef Hidden_type Algebraic_kernel_1; 

/*! 
represents coordinates of points. 
Typedef from `Algebraic_kernel_1::Algebraic_real_1` 
*/ 
typedef Hidden_type Algebraic_real_1; 

/*! 
Typedef from `Algebraic_kernel_1::Bound` 
*/ 
typedef Hidden_type Bound; 

/// @} 

/// \name Accessing functor objects 
/// @{

/*! 

*/ 
Construct_curve_2 construct_curve_2_object() const; 

/*! 

*/ 
Construct_point_2 construct_point_2_object() const; 

/*! 

*/ 
Construct_x_monotone_segment_2 construct_x_monotone_segment_2_object() const; 

/// @}

}; /* end Arr_algebraic_segment_traits_2 */
} /* end namespace CGAL */
