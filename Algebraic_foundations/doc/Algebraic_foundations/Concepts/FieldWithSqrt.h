
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

A model of `FieldWithSqrt` is a `Field` that has operations to take square roots.

Moreover, `CGAL::Algebraic_structure_traits< FieldWithSqrt >` is a model of `AlgebraicStructureTraits` providing:

- \link AlgebraicStructureTraits::Algebraic_category `CGAL::Algebraic_structure_traits< FieldWithSqrt >::Algebraic_category` \endlink derived from `CGAL::Field_with_sqrt_tag`
- \link AlgebraicStructureTraits::Sqrt `CGAL::Algebraic_structure_traits< FieldWithSqrt >::Sqrt` \endlink which is a model of `AlgebraicStructureTraits_::Sqrt`

\cgalRefines `Field`

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

class FieldWithSqrt {
public:

}; /* end FieldWithSqrt */

