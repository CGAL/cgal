
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableFunctor` provides evaluation of a
`PolynomialTraits_d::Polynomial_d` interpreted as a homogeneous polynomial
<B>in one variable</B>.

For instance the polynomial \f$ p = 5x^2y^3 + y\f$ is interpreted as the homogeneous polynomial
\f$ p[x](u,v) = 5x^2u^3 + uv^2\f$ and evaluated as such.

\cgalRefines `AdaptableFunctor`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::EvaluateHomogeneous {
public:

/// \name Types
/// @{

/*!

*/
typedef PolynomialTraits_d::Coefficient_type result_type;

/// @}

/// \name Operations
/// @{

/*!

Returns \f$ p(u,v)\f$, with respect to the outermost variable.

*/
result_type operator()(PolynomialTraits_d::Polynomial_d p,
PolynomialTraits_d::Coefficient_type u,
PolynomialTraits_d::Coefficient_type v);

/// @}

}; /* end PolynomialTraits_d::EvaluateHomogeneous */

