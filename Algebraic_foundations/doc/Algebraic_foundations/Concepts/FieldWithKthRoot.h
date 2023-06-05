
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

A model of `FieldWithKthRoot` is a `FieldWithSqrt` that has operations to take k-th roots.

Moreover, `CGAL::Algebraic_structure_traits< FieldWithKthRoot >` is a model of `AlgebraicStructureTraits` providing:

- \link AlgebraicStructureTraits::Algebraic_category `CGAL::Algebraic_structure_traits< FieldWithKthRoot >::Algebraic_category` \endlink derived from `CGAL::Field_with_kth_root_tag`
- \link AlgebraicStructureTraits::Kth_root `CGAL::Algebraic_structure_traits< FieldWithKthRoot >::Kth_root` \endlink which is a model of `AlgebraicStructureTraits_::KthRoot`

\cgalRefines{FieldWithSqrt}

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

class FieldWithKthRoot {

}; /* end FieldWithKthRoot */

