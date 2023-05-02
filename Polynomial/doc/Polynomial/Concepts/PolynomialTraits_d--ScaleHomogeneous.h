
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

Given a numerator \f$ a\f$ and a denominator \f$ b\f$ this `AdaptableFunctor`
scales a `PolynomialTraits_d::Polynomial_d` \f$ p\f$ with respect to one variable,
that is, it computes \f$ b^{degree(p)}\cdot p(a/b\cdot x)\f$.

Note that this functor operates on the polynomial in the univariate view, that is,
the polynomial is considered as a univariate homogeneous polynomial in one specific variable.

\cgalRefines{AdaptableFunctor,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::ScaleHomogeneous {
public:

/// \name Types
/// @{

/*!

*/
typedef PolynomialTraits_d::Polynomial_d result_type;

/// @}

/// \name Operations
/// @{

/*!
Returns \f$ b^{degree}\cdot p(a/b\cdot x)\f$,
with respect to the outermost variable.
*/
result_type operator()(PolynomialTraits_d::Polynomial_d p,
PolynomialTraits_d::Innermost_coefficient_type a,
PolynomialTraits_d::Innermost_coefficient_type b);

/*!
Same as first operator but for variable \f$ x_i\f$.
\pre \f$ 0 \leq i < d\f$.

*/
result_type operator()(PolynomialTraits_d::Polynomial_d p,
PolynomialTraits_d::Innermost_coefficient_type a,
PolynomialTraits_d::Innermost_coefficient_type b,
int i);

/// @}

}; /* end PolynomialTraits_d::ScaleHomogeneous */

