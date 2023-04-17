
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableFunctor` computes the <I>pseudo division</I>
of two polynomials \f$ f\f$ and \f$ g\f$.

Given \f$ f\f$ and \f$ g \neq 0\f$ this functor computes quotient \f$ q\f$ and
remainder \f$ r\f$ such that \f$ D \cdot f = g \cdot q + r\f$ and \f$ degree(r) < degree(g)\f$,
where \f$ D = leading\_coefficient(g)^{max(0, degree(f)-degree(g)+1)}\f$

This functor is useful if the regular division is not available,
which is the case if `PolynomialTraits_d::Coefficient_type` is not a `Field`.
Hence in general it is not possible to invert the leading coefficient of \f$ g\f$.
Instead \f$ f\f$ is extended by \f$ D\f$ allowing integral divisions in the internal
computation.

\cgalRefines{AdaptableFunctor,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::PseudoDivision`
\sa `PolynomialTraits_d::PseudoDivisionRemainder`
\sa `PolynomialTraits_d::PseudoDivisionQuotient`

*/

class PolynomialTraits_d::PseudoDivision {
public:

/// \name Types
/// @{

/*!

*/
typedef void result_type;

/// @}

/// \name Operations
/// @{

/*!

Computes the pseudo division with respect to the outermost variable
\f$ x_{d-1}\f$.

*/
result_type operator()(PolynomialTraits_d::Polynomial_d f,
PolynomialTraits_d::Polynomial_d g,
PolynomialTraits_d::Polynomial_d & q,
PolynomialTraits_d::Polynomial_d & r,
PolynomialTraits_d::Coefficient_type & D);

/// @}

}; /* end PolynomialTraits_d::PseudoDivision */

