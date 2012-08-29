/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  A model of `UniqueFactorizationDomain` is an `IntegralDomain` with the 
///  additional property 
///  that the ring it represents is a unique factorization domain 
///  (a.k.a. UFD or factorial ring), meaning that every non-zero non-unit 
///  element has a factorization into irreducible elements that is unique 
///  up to order and up to multiplication by invertible elements (units). 
///  (An irreducible element is a non-unit ring element that cannot be factored 
///  further into two non-unit elements. In a UFD, the irreducible elements 
///  are precisely the prime elements.)
///  In a UFD, any two elements, not both zero, possess a greatest common 
///  divisor (gcd).
///  Moreover, `CGAL::Algebraic_structure_traits< UniqueFactorizationDomain >` 
///  is a model of `AlgebraicStructureTraits` providing: 
///  - `CGAL::Algebraic_structure_traits::Algebraic_type` 
///  derived from `Unique_factorization_domain_tag` 
///  - `CGAL::Algebraic_structure_traits::Gcd`
///  \refines ::IntegralDomain
///  \sa `IntegralDomainWithoutDivision`
///  \sa `IntegralDomain`
///  \sa `UniqueFactorizationDomain`
///  \sa `EuclideanRing`
///  \sa `Field`
///  \sa `FieldWithSqrt`
///  \sa `FieldWithKthRoot`
///  \sa `FieldWithRootOf`
///  \sa `AlgebraicStructureTraits`
class UniqueFactorizationDomain {
public:

}; /* concept UniqueFactorizationDomain */
/// @}
/// @} 

                   
  

