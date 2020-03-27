
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the total degree
of a `PolynomialTraits_d::Polynomial_d`.

Given a (multivariate) monomial the sum of all appearing exponents
is the total degree of this monomial.
The total degree of a polynomial \f$ p\f$ is the maximum of the total degrees
of all appearing (multivariate) monomials in \f$ p\f$.

For instance the total degree of \f$ p = x_0^2x_1^3+x_1^4\f$ is \f$ 5\f$.

The total degree of the zero polynomial is set to \f$ 0\f$.
From the mathematical point of view this should
be \f$ -\infty\f$, but this would imply an inconvenient return type.

\cgalRefines `AdaptableUnaryFunction`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::Degree`
\sa `PolynomialTraits_d::DegreeVector`

*/

class PolynomialTraits_d::TotalDegree {
public:

/// \name Types
/// @{

/*!

*/
typedef int result_type;

/*!

*/
typedef PolynomialTraits_d::Polynomial_d argument_type;

/// @}

/// \name Operations
/// @{

/*!
Computes the total degree of \f$ p\f$.
*/
result_type operator()(argument_type p);

/// @}

}; /* end PolynomialTraits_d::TotalDegree */

