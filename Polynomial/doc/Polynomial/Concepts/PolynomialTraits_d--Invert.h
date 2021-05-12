
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` inverts one variable in a given
`PolynomialTraits_d::Polynomial_d`, that is, for a given polynomial
\f$ p\f$ it computes \f$ x^{degree(p)}p(1/x)\f$.

Note that this functor operates on the polynomial in the univariate view, that is,
the polynomial is considered as a univariate polynomial in one specific variable.

This functor is provided for efficiency reasons, since this operation just inverts the
order of the coefficients with respect to the specified variable.

\cgalRefines `AdaptableUnaryFunction`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::Invert {
public:

/// \name Types
/// @{

/*!

*/
typedef PolynomialTraits_d::Polynomial_d result_type;

/*!

*/
typedef PolynomialTraits_d::Polynomial_d argument_type;

/// @}

/// \name Operations
/// @{

/*!

Returns \f$ x^{degree(p)}p(1/x)\f$,
where x refers to the outermost variable \f$ x_{d-1}\f$.

*/
result_type operator()(argument_type p);

/*!

Return \f$ x^{degree(p,i)}p(1/x)\f$,
where x refers to the variable \f$ x_{i}\f$.
\pre \f$ 0 \leq i < d\f$.
*/
result_type operator()(argument_type p, int i);

/// @}

}; /* end PolynomialTraits_d::Invert */

