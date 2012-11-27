
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

A model of `FieldWithRootOf` is a `FieldWithKthRoot` with the possibility to 
construct it as the root of a univariate polynomial. 

Moreover, `CGAL::Algebraic_structure_traits< FieldWithRootOf >` is a model of `AlgebraicStructureTraits` providing: 

- `CGAL::Algebraic_structure_traits< FieldWithRootOf >::Algebraic_type` derived from `CGAL::Field_with_kth_root_tag` 
- `CGAL::Algebraic_structure_traits< FieldWithRootOf >::Root_of`  which is a model of `AlgebraicStructureTraits::RootOf`

\cgalRefines `FieldWithKthRoot` 

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

class FieldWithRootOf {
public:

/// @}

}; /* end FieldWithRootOf */

