
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

For a given `PolynomialTraits_d::Polynomial_d` \f$ p\f$
this `AdaptableUnaryFunction` returns the degree vector, that is,
it returns the exponent vector of the monomial of highest order in \f$ p\f$,
where the monomial order is the lexicographic order giving outer
variables a higher priority. In particular, this is the monomial
that belongs to the innermost leading coefficient of \f$ p\f$.

\cgalRefines `AdaptableUnaryFunction`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::Degree`
\sa `PolynomialTraits_d::TotalDegree`
\sa `PolynomialTraits_d::InnermostLeadingCoefficient`

*/

class PolynomialTraits_d::DegreeVector {
public:

/// \name Types
/// @{

/*!

*/
  typedef CGAL::Exponent_vector result_type;

/*!

*/
typedef PolynomialTraits_d::Polynomial_d argument_type;

/// @}

/// \name Operations
/// @{

/*!
Returns the degree vector.
*/
result_type operator()(argument_type p);

/// @}

}; /* end PolynomialTraits_d::DegreeVector */

