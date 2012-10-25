
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalconcept

A model of `FieldWithSqrt` is a `Field` that has operations to take square roots. 

Moreover, `CGAL::Algebraic_structure_traits< FieldWithSqrt >` is a model of `AlgebraicStructureTraits` providing: 

- `CGAL::Algebraic_structure_traits< FieldWithSqrt >::Algebraic_type` derived from `CGAL::Field_with_sqrt_tag` 
- `CGAL::Algebraic_structure_traits< FieldWithSqrt >::Sqrt`  which is a model of `AlgebraicStructureTraits::Sqrt` 

\refines `Field` 

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

class FieldWithSqrt {
public:

/// @}

}; /* end FieldWithSqrt */

