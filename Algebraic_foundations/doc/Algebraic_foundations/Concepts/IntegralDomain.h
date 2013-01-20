
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`IntegralDomain` refines `IntegralDomainWithoutDivision` by 
providing an integral division. 

<B>Note:</B> The concept does not require the operator / for this operation. 
We intend to reserve the operator syntax for use with a `Field`. 

Moreover, `CGAL::Algebraic_structure_traits< IntegralDomain >` is a model of 
`AlgebraicStructureTraits` providing: 

- \link AlgebraicStructureTraits::Algebraic_category `CGAL::Algebraic_structure_traits< IntegralDomain >::Algebraic_category` \endlink derived from `CGAL::Integral_domain_tag` 
- \link AlgebraicStructureTraits::Integral_division `CGAL::Algebraic_structure_traits< IntegralDomain >::Integral_division` \endlink  which is a model of `AlgebraicStructureTraits::IntegralDivision`
- `CGAL::Algebraic_structure_traits< IntegralDomain >::Divides`  which is a model of `AlgebraicStructureTraits::Divides`

\cgalRefines `IntegralDomainWithoutDivision` 

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

class IntegralDomain {
public:

/// @}

}; /* end IntegralDomain */

