
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

This is the most basic concept for algebraic structures considered within \cgal.

A model `IntegralDomainWithoutDivision` represents an integral domain,
i.e.\ commutative ring with 0, 1, +, * and unity free of zero divisors.

<B>Note:</B> A model is not required to offer the always well defined integral division.

It refines `Assignable`, `CopyConstructible`, `DefaultConstructible`
and `FromIntConstructible`.

It refines `EqualityComparable`, where equality is defined w.r.t.\ the ring element being represented.

The operators unary and binary plus +, unary and binary minus -,
multiplication * and their compound forms +=, -=, *= are required and
implement the respective ring operations.

Moreover, `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >` is a model of
`AlgebraicStructureTraits` providing:

- \link AlgebraicStructureTraits::Algebraic_category `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >::Algebraic_category` \endlink derived from `CGAL::Integral_domain_without_division_tag`
- \link AlgebraicStructureTraits::Is_zero `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >::Is_zero` \endlink  which is a model of `AlgebraicStructureTraits_::IsZero`
- \link AlgebraicStructureTraits::Is_one `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >::Is_one` \endlink  which is a model of `AlgebraicStructureTraits_::IsOne`
- \link AlgebraicStructureTraits::Square `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >::Square` \endlink  which is a model of `AlgebraicStructureTraits_::Square`
- \link AlgebraicStructureTraits::Simplify `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >::Simplify` \endlink which is a model of `AlgebraicStructureTraits_::Simplify`
- \link AlgebraicStructureTraits::Unit_part `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >::Unit_part` \endlink  which is a model of `AlgebraicStructureTraits_::UnitPart`

\cgalRefines `Assignable`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`
\cgalRefines `EqualityComparable`
\cgalRefines `FromIntConstructible`

\sa `IntegralDomainWithoutDivision`
\sa `IntegralDomain`
\sa `UniqueFactorizationDomain`
\sa `EuclideanRing`
\sa `Field`
\sa `FieldWithSqrt`
\sa `FieldWithKthRoot`
\sa `FieldWithRootOf`
\sa `AlgebraicStructureTraits`

*/

class IntegralDomainWithoutDivision {
public:

/// \name Operations
/// @{

/*!
unary plus
*/
IntegralDomainWithoutDivision
operator+(const IntegralDomainWithoutDivision &a);

/*!
unary minus
*/
IntegralDomainWithoutDivision
operator-(const IntegralDomainWithoutDivision &a);

/*!

*/
IntegralDomainWithoutDivision
operator+(const IntegralDomainWithoutDivision &a,
const IntegralDomainWithoutDivision &b);

/*!

*/
IntegralDomainWithoutDivision
operator-(const IntegralDomainWithoutDivision &a,
const IntegralDomainWithoutDivision &b);

/*!

*/
IntegralDomainWithoutDivision
operator*(const IntegralDomainWithoutDivision &a,
const IntegralDomainWithoutDivision &b);


/*!

*/
IntegralDomainWithoutDivision
operator+=(const IntegralDomainWithoutDivision &b);

/*!

*/
IntegralDomainWithoutDivision
operator-=(const IntegralDomainWithoutDivision &b);

/*!

*/
IntegralDomainWithoutDivision
operator*=(const IntegralDomainWithoutDivision &b);

/*!
The `result_type` is convertible to `bool`.
*/
result_type
operator==(const IntegralDomainWithoutDivision &a, const IntegralDomainWithoutDivision &b);

/*!
The `result_type` is convertible to `bool`.
*/
result_type
operator!=(const IntegralDomainWithoutDivision &a, const IntegralDomainWithoutDivision &b);

/// @}

}; /* end IntegralDomainWithoutDivision */

