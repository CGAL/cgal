
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableBinaryFunction` compares two polynomials, with respect to the lexicographic
order with preference to the outermost variable.

This functor is well defined if `PolynomialTraits_d::Innermost_coefficient_type` is
`RealEmbeddable`.

\cgalRefines{AdaptableBinaryFunction,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::Compare {
public:

/// \name Types
/// @{

/*!

*/
typedef CGAL::Comparison_result result_type;

/*!

*/
typedef PolynomialTraits_d::Polynomial_d first_argument_type;

/*!

*/
typedef PolynomialTraits_d::Polynomial_d second_argument_type;

/// @}

/// \name Operations
/// @{

/*!
Compares two polynomials.
*/
result_type operator()(first_argument_type f,
second_argument_type g);

/// @}

}; /* end PolynomialTraits_d::Compare */

