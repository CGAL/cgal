
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

For the given `PolynomialTraits_d::Polynomial_d` this
`AdaptableBinaryFunction` returns the coefficient of
the (multivariate) monomial specified by the given `CGAL::Exponent_vector`.

\cgalRefines{AdaptableBinaryFunction,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::GetInnermostCoefficient {
public:

/// \name Types
/// @{

/*!

*/
typedef PolynomialTraits_d::Innermost_coefficient_type result_type;

/*!

*/
typedef PolynomialTraits_d::Polynomial_d first_argument_type ;

/*!

*/
  typedef CGAL::Exponent_vector second_argument_type;

/// @}

/// \name Operations
/// @{

/*!

For given polynomial \f$ p\f$ this operator returns the innermost coefficient of the
monomial corresponding to the given `CGAL::Exponent_vector` \f$ v\f$.
*/
result_type operator()( first_argument_type p,
second_argument_type v);

/// @}

}; /* end PolynomialTraits_d::GetInnermostCoefficient */

