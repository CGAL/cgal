
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

A model of `Field` is an `IntegralDomain` in which every non-zero element
has a multiplicative inverse.
Thus, one can divide by any non-zero element.
Hence division is defined for any divisor != 0.
For a Field, we require this division operation to be available through
operators / and /=.

Moreover, `CGAL::Algebraic_structure_traits< Field >` is a model of
`AlgebraicStructureTraits` providing:

- \link AlgebraicStructureTraits::Algebraic_category `CGAL::Algebraic_structure_traits< Field >::Algebraic_category` \endlink derived from `CGAL::Field_tag`
- \link AlgebraicStructureTraits::Inverse `CGAL::Algebraic_structure_traits< FieldWithSqrt >::Inverse` \endlink  which is a model of `AlgebraicStructureTraits_::Inverse`

\cgalRefines{IntegralDomain}

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

class Field {
public:

/// \name Operations
/// @{

/*!

*/
Field operator/(const Field &a, const Field &b);


/*!

*/
Field operator/=(const Field &b);

/// @}

}; /* end Field */

