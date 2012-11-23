
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

A model of `FieldWithKthRoot` is a `FieldWithSqrt` that has operations to take k-th roots. 

Moreover, `CGAL::Algebraic_structure_traits< FieldWithKthRoot >` is a model of `AlgebraicStructureTraits` providing: 

- `CGAL::Algebraic_structure_traits< FieldWithKthRoot >::Algebraic_type` derived from `Field_with_kth_root_tag` 

- `CGAL::Algebraic_structure_traits< FieldWithKthRoot >::Kth_root` which is a model of `AlgebraicStructureTraits::KthRoot`

\cgalRefines `FieldWithSqrt` 

\sa ::IntegralDomainWithoutDivision 
\sa ::IntegralDomain 
\sa ::UniqueFactorizationDomain 
\sa ::EuclideanRing 
\sa ::Field 
\sa ::FieldWithSqrt 
\sa ::FieldWithKthRoot 
\sa ::FieldWithRootOf 
\sa ::AlgebraicStructureTraits 

*/

class FieldWithKthRoot {

}; /* end FieldWithKthRoot */

