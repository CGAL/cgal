
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

For a given polynomial \f$ p\f$ this `AdaptableUnaryFunction` computes the
unique representative of the set
\f[ {\cal P} := \{ q\ |\ \lambda * q = p\ for\ some\ \lambda \in R \}, \f]
where \f$ R\f$ is the base of the polynomial ring.

In case `PolynomialTraits::Innermost_coefficient_type` is a model of
`Field`, the computed polynomial is the <I>monic</I> polynomial in
\f$ \cal P\f$, that is, the innermost leading coefficient equals one.

In case `PolynomialTraits::Innermost_coefficient_type` is a model
of `UniqueFactorizationDomain`, the computed polynomial is the one with
a multivariate content of one.

For all other cases the notion of uniqueness is up to the concrete model.

Note that the computed polynomial has the same zero set as the given one.

\cgalRefines{AdaptableUnaryFunction,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::Canonicalize {
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

Returns the canonical representative of \f$ p\f$.
*/
result_type operator()(first_argument_type p);

/// @}

}; /* end PolynomialTraits_d::Canonicalize */

