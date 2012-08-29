/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  ::IntegralDomain refines ::IntegralDomainWithoutDivision by 
///  providing an integral division.
///  <B>Note:</B> The concept does not require the operator / for this operation. 
///  We intend to reserve the operator syntax for use with a `Field`.
///  
///  
///  Moreover, `CGAL::Algebraic_structure_traits< IntegralDomain >` is a model of 
///  `AlgebraicStructureTraits` providing:
///  - `CGAL::Algebraic_structure_traits::Algebraic_type` derived from `Integral_domain_tag` 
///  - `CGAL::Algebraic_structure_traits::Integral_division`
///  - `CGAL::Algebraic_structure_traits::Divides`
///
///  \refines ::IntegralDomainWithoutDivision
///  \sa `IntegralDomainWithoutDivision`
///  \sa `IntegralDomain`
///  \sa `UniqueFactorizationDomain`
///  \sa `EuclideanRing`
///  \sa `Field`
///  \sa `FieldWithSqrt`
///  \sa `FieldWithKthRoot`
///  \sa `FieldWithRootOf`
///  \sa `AlgebraicStructureTraits`
class IntegralDomain {
public:

}; /* concept IntegralDomain */
/// @}
/// @} 

                   
  

