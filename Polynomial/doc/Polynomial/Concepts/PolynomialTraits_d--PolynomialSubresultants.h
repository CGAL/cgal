
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

Computes the polynomial subresultant of two polynomials \f$ p\f$ and \f$ q\f$ of
type `PolynomialTraits_d::Polynomial_d` with respect to outermost variable.
Let
\f$ p=\ccSum{i=0,\ldots,n}{} p_i x^i\f$ and
\f$ q=\ccSum{i=0,\ldots,m}{} q_i x^i\f$, where \f$ x\f$
is the outermost variable.
The \f$ i\f$-th subresultant (with \f$ i=0,\ldots,\min\{n,m\}\f$) is defined by

\f[
\mathrm{Sres}_{i}(p,q) = \det
\begin{pmatrix}
    p_{n} & \dots & & \dots & p_{2i-m+2} & x^{m-i-1}p  \\
    & \ddots & & & \vdots & \vdots\\
    & & p_{n} & \dots & p_{i+1} & p \\
    q_{m} & \dots & & \dots & q_{2i-n+2} & x^{n-i-1}q  \\
    & \ddots & & & \vdots & \vdots\\
    & & q_{m} & \dots & q_{i+1} & q
\end{pmatrix}
\f]

where \f$ p_i\f$ and \f$ q_i\f$ are set to zero if \f$ i<0\f$.
In the case that \f$ n=m\f$, \f$ \mathrm{Sres_n}\f$ is set to \f$ q\f$.

The result is written in an output range, starting with the \f$ 0\f$-th subresultant
\f$ \mathrm{Sres}_0(p,q)\f$
(aka as the resultant of \f$ p\f$ and \f$ q\f$).

\note This functor is optional.

\cgalRefines{AdaptableBinaryFunction,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::Resultant`
\sa `PolynomialTraits_d::PrincipalSubresultants`
\sa `PolynomialTraits_d::PolynomialSubresultantsWithCofactors`
\sa `PolynomialTraits_d::SturmHabichtSequence`

*/

class PolynomialTraits_d::PolynomialSubresultants {
public:

/// \name Operations
/// @{

/*!
computes the polynomial subresultants of \f$ p\f$ and \f$ q\f$,
with respect to the outermost variable. Each element is of type
`PolynomialTraits_d::Polynomial_d`.
*/
template<typename OutputIterator>
OutputIterator operator()(Polynomial_d p,
Polynomial_d q,
OutputIterator out);

/*!
computes the polynomial subresultants of \f$ p\f$ and \f$ q\f$,
with respect to the variable \f$ x_i\f$.
*/
template<typename OutputIterator>
OutputIterator operator()(Polynomial_d p,
Polynomial_d q,
OutputIterator out,
int i);

/// @}

}; /* end PolynomialTraits_d::PolynomialSubresultants */

