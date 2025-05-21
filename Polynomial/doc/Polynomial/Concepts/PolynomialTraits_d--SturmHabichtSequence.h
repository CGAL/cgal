
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

Computes the Sturm-Habicht sequence
(aka the signed subresultant sequence)
of a polynomial \f$ f\f$ of type
`PolynomialTraits_d::Polynomial_d` with respect to a certain variable \f$ x_i\f$.
The Sturm-Habicht sequence is similar to the polynomial subresultant sequence
of \f$ f\f$ and its derivative \f$ f':=\frac{\partial f}{\partial x_i}\f$
with respect to \f$ x_i\f$. The implementation is based on the following definition:

Let \f$ n:=\deg f\f$ and \f$ \delta_k:=(-1)^{k(k+1)/2}\f$.
For \f$ k\in\{0,\ldots,n\}\f$, the <I>\f$ k\f$-th Sturm-Habicht polynomial</I>
of \f$ f\f$ is defined as:

\f[
\mathrm{Stha}_{k}(f) = \left \{
\begin{array}{cll}
f & \text{if} & k = n \\
f' & \text{if} & k = n - 1 \\
\delta_{n - k - 1}\mathrm{Sres}_{k}(f,f') & \text{if} & 0 \leq k \leq n - 2
\end{array}
\right .
\f]

where \f$ \mathrm{Sres}_k(f,f')\f$ is defined
as in the concept `PolynomialTraits_d::PolynomialSubresultants`.

The result is written in an output range,
starting with the \f$ 0\f$-th Sturm-Habicht polynomial (which is equal to
the discriminant of \f$ f\f$ up to a multiple of the leading coefficient).

\note This functor is optional.

\cgalRefines{AdaptableBinaryFunction,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::Resultant`
\sa `PolynomialTraits_d::PrincipalSturmHabichtSequence`
\sa `PolynomialTraits_d::SturmHabichtSequenceWithCofactors`
\sa `PolynomialTraits_d::PolynomialSubresultants`

*/

class PolynomialTraits_d::SturmHabichtSequence {
public:

/// \name Operations
/// @{

/*!
computes the Sturm-Habicht sequence of \f$ f\f$,
with respect to the outermost variable. Each element is of type
`PolynomialTraits_d::Polynomial_d`.
*/
template<typename OutputIterator>
OutputIterator operator()(Polynomial_d f,
OutputIterator out);

/*!
computes the Sturm-Habicht sequence of \f$ f\f$
with respect to the variable \f$ x_i\f$.
*/
template<typename OutputIterator>
OutputIterator operator()(Polynomial_d f,
OutputIterator out,
int i);

/// @}

}; /* end PolynomialTraits_d::SturmHabichtSequence */

