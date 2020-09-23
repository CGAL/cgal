
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

Computes the principal subresultant of two polynomials \f$ p\f$ and \f$ q\f$ of
type `PolynomialTraits_d::Coefficient_type`
with respect to the outermost variable.
The \f$ i\f$-th principal subresultant, \f$ \mathrm{sres}_i(p,q)\f$,
is defined as the coefficient at \f$ t^i\f$ of the \f$ i\f$-th polynomial
subresultant \f$ \mathrm{Sres}_i(p,q)\f$. Thus, it is either the leading
coefficient of \f$ \mathrm{Sres}_i\f$, or zero in the case where its degree is
below \f$ i\f$.

The result is written in an output range, starting with the \f$ 0\f$-th
principal subresultant \f$ \mathrm{sres}_0(p,q)\f$
,aka as the resultant of \f$ p\f$ and \f$ q\f$.
(Note that \f$ \mathrm{sres}_0(p,q)=\mathrm{Sres}_0(p,q)\f$ by definition)

\note This functor is optional.

\cgalRefines `AdaptableBinaryFunction`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::Resultant`
\sa `PolynomialTraits_d::PolynomialSubresultants`
\sa `PolynomialTraits_d::PrincipalSturmHabichtSequence`

*/

class PolynomialTraits_d::PrincipalSubresultants {
public:

/// \name Operations
/// @{

/*!
computes the principal subresultants of \f$ p\f$ and \f$ q\f$,
with respect to the outermost variable. Each element is of type
`PolynomialTraits_d::Coefficient_type`.
*/
template<typename OutputIterator>
OutputIterator operator()(Polynomial_d p,
Polynomial_d q,
OutputIterator out);

/*!
computes the principal subresultants of \f$ p\f$ and \f$ q\f$,
with respect to the variable \f$ x_i\f$.
*/
template<typename OutputIterator>
OutputIterator operator()(Polynomial_d p,
Polynomial_d q,
OutputIterator out,
int i);

/// @}

}; /* end PolynomialTraits_d::PrincipalSubresultants */

