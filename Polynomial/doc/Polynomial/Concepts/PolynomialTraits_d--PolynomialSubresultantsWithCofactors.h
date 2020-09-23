
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

Computes the polynomial subresultant of two polynomials \f$ p\f$ and \f$ q\f$ of degree
\f$ n\f$ and \f$ m\f$, respectively,
as defined in the documentation of `PolynomialTraits_d::PolynomialSubresultants`.
Moreover, for \f$ \mathrm{Sres}_i(p,q)\f$, polynomials \f$ u_i\f$ and \f$ v_i\f$
with \f$ \deg u_i\leq m-i-1\f$ and \f$ \deg v_i\leq n-i-1\f$ are computed
such that \f$ \mathrm{Sres}_i(p,q)=u_i p + v_i q\f$. \f$ u_i\f$ and \f$ v_i\f$ are called
the <I>cofactors</I> of \f$ \mathrm{Sres}_i(p,q)\f$.

The result is written in three output ranges, each of length \f$ \min\{n,m\}+1\f$,
starting with the \f$ 0\f$-th subresultant and the corresponding cofactors.

\note This functor is optional.

\cgalRefines `AdaptableBinaryFunction`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::Resultant`
\sa `PolynomialTraits_d::PolynomialSubresultants`
\sa `PolynomialTraits_d::PrincipalSubresultants`
\sa `PolynomialTraits_d::SturmHabichtSequenceWithCofactors`

*/

class PolynomialTraits_d::PolynomialSubresultantsWithCofactors {
public:

/// \name Operations
/// @{

/*!
computes the subresultants of \f$ p\f$ and \f$ q\f$, and the cofactors,
with respect to the outermost variable. Each element is of type
`PolynomialTraits_d::Polynomial_d`.
*/
template< typename OutputIterator1,
typename OutputIterator2,
typename OutputIterator3 >
OutputIterator1 operator()(Polynomial_d p,
Polynomial_d q,
OutputIterator1 sres,
OutputIterator2 co_p,
OutputIterator3 co_q);

/*!
computes the subresultants of \f$ p\f$ and \f$ q\f$, and the cofactors,
with respect to \f$ x_i\f$. Each element is of type
`PolynomialTraits_d::Polynomial_d`.
*/
template< typename OutputIterator1,
typename OutputIterator2,
typename OutputIterator3 >
OutputIterator1 operator()(Polynomial_d p,
Polynomial_d q,
OutputIterator1 sres,
OutputIterator2 co_p,
OutputIterator3 co_q,
int i);

/// @}

}; /* end PolynomialTraits_d::PolynomialSubresultantsWithCofactors */

