
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableFunctor` moves a variable at position \f$ i\f$ to a
new position \f$ j\f$. The relative order of the other variables is preserved,
that is, the variables between \f$ x_i\f$ and \f$ x_j\f$ (including \f$ x_j\f$) are moved by
one position while \f$ x_i\f$ is moved to the former position of \f$ x_j\f$.

This function may be used to make a certain variable the outer most variable.

\cgalRefines `AdaptableFunctor`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::Move {
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

This function moves the variable at position \f$ i\f$ to its new position \f$ j\f$ and returns
the new polynomial. The relative order of the other variables is preserved.
\pre \f$ 0 \leq i < d\f$.
\pre \f$ 0 \leq j < d\f$.

*/
result_type operator()(PolynomialTraits_d::Polynomial_d,
int i, int j);

/// @}

}; /* end PolynomialTraits_d::Move */

