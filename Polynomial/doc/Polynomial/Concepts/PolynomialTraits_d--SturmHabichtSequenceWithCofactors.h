
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

Computes the Sturm-Habicht polynomials of a polynomial \f$ f\f$ of degree \f$ n\f$,
as defined in the documentation of `PolynomialTraits_d::SturmHabichtSequence`.
Moreover, for \f$ \mathrm{Stha}_i(f)\f$, polynomials \f$ u_i\f$ and \f$ v_i\f$
with \f$ \deg u_i\leq n-i-2\f$ and \f$ \deg v_i\leq n-i-1\f$ are computed
such that \f$ \mathrm{Sres}_i(p,q)=u_i f + v_i f'\f$. \f$ u_i\f$ and \f$ v_i\f$ are called
the <I>cofactors</I> of \f$ \mathrm{Stha}_i(f)\f$.

The result is written in three output ranges, each of length \f$ \min\{n,m\}+1\f$,
starting with the \f$ 0\f$-th Sturm-Habicht polynomial \f$ \mathrm{Stha_0(f)}\f$
and the corresponding cofactors.

\note This functor is optional.

\cgalRefines `AdaptableBinaryFunction`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::Resultant`
\sa `PolynomialTraits_d::SturmHabichtSequence`
\sa `PolynomialTraits_d::PrincipalSturmHabichtSequence`
\sa `PolynomialTraits_d::PolynomialSubresultantsWithCofactors`

*/

class PolynomialTraits_d::SturmHabichtSequenceWithCofactors {
public:

/// \name Operations
/// @{

/*!
computes the Sturm-Habicht sequence of \f$ f\f$, and the cofactors,
with respect to the outermost variable. Each element is of type
`PolynomialTraits_d::Polynomial_d`.
*/
template<typename OutputIterator1,
typename OutputIterator2,
typename OutputIterator3>
OutputIterator1 operator()(Polynomial_d f,
OutputIterator1 stha,
OutputIterator2 co_f,
OutputIterator3 co_fx);

/*!
computes the Sturm-Habicht sequence of \f$ f\f$, and the cofactors,
with respect to \f$ x_i\f$. Each element is of type
`PolynomialTraits_d::Polynomial_d`.
*/
template< typename OutputIterator1,
typename OutputIterator2,
typename OutputIterator3 >
OutputIterator1 operator()(Polynomial_d f,
OutputIterator1 stha,
OutputIterator2 co_f,
OutputIterator3 co_fx,
int i);

/// @}

}; /* end PolynomialTraits_d::SturmHabichtSequenceWithCofactors */

