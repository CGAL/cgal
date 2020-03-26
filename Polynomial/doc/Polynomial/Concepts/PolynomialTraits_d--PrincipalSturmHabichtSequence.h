
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

Computes the principal leading coefficients of the Sturm-Habicht sequence
of a polynomials \f$ f\f$ of type `PolynomialTraits_d::Polynomial_d`
with respect a certain variable \f$ x_i\f$.
This means that for the \f$ j\f$-th Sturm-Habicht polynomial, this methods returns
the coefficient of \f$ x_i^j\f$.

Note that the degree of the \f$ j\f$-th Sturm-Habicht polynomial is at most \f$ j\f$,
but the principal coefficient might be zero, thus, this functor does not
necessarily give the leading coefficient of the Sturm-Habicht polynomials.

In case that `PolynomialTraits_d::Coefficient_type` is `RealEmbeddable`, the function `CGAL::number_of_real_roots` can be used
on the resulting sequence to count the number of distinct real roots of
the polynomial \f$ f\f$.

\note This functor is optional.

\cgalRefines `AdaptableBinaryFunction`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `CGAL::number_of_real_roots()`
\sa `PolynomialTraits_d::Resultant`
\sa `PolynomialTraits_d::SturmHabichtSequence`
\sa `PolynomialTraits_d::PrincipalSubresultants`

*/

class PolynomialTraits_d::PrincipalSturmHabichtSequence {
public:

/// \name Operations
/// @{

/*!
computes the principal coefficients of the
Sturm-Habicht sequence of \f$ f\f$,
with respect to the outermost variable. Each element is of type
`PolynomialTraits_d::Coefficient_type`.
*/
template<typename OutputIterator>
OutputIterator operator()(Polynomial_d f,
OutputIterator out);

/*!
computes the principal coefficients
of the Sturm-Habicht sequence of \f$ f\f$
with respect to the variable \f$ x_i\f$.
*/
template<typename OutputIterator>
OutputIterator operator()(Polynomial_d f,
OutputIterator out,
int i);

/// @}

}; /* end PolynomialTraits_d::PrincipalSturmHabichtSequence */

