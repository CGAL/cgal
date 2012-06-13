/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  A model of `EuclideanRing` represents an euclidean ring (or Euclidean domain). 
///  It is an `UniqueFactorizationDomain` that affords a suitable notion of minimality of remainders 
///  such that given \f$x\f$ and \f$y  \neq 0\f$ we obtain an (almost) unique solution to 
///  \f$ x = qy + r \f$ by demanding that a solution \f$(q,r)\f$ is chosen to minimize \f$r\f$. 
///  In particular, \f$r\f$ is chosen to be \f$0\f$ if possible.
///  Moreover, `CGAL::Algebraic_structure_traits< EuclideanRing >` is a model of 
///  `AlgebraicStructureTraits` providing:
///  - `CGAL::Algebraic_structure_traits::Algebraic_type` derived from `Unique_factorization_domain_tag` 
///  - `CGAL::Algebraic_structure_traits::Mod` 
///  - `CGAL::Algebraic_structure_traits::Div` 
///  - `CGAL::Algebraic_structure_traits::Div_mod`
///  Remarks 
///  -------------- 
///  
///  The most prominent example of a Euclidean ring are the integers. 
///  Whenever both \f$x\f$ and \f$y\f$ are positive, then it is conventional to choose 
///  the smallest positive remainder \f$r\f$.
///  
///  
///  
///  
///  \refines ::UniqueFactorizationDomain
///  \sa `IntegralDomainWithoutDivision`
///  \sa `IntegralDomain`
///  \sa `UniqueFactorizationDomain`
///  \sa `EuclideanRing`
///  \sa `Field`
///  \sa `FieldWithSqrt`
///  \sa `FieldWithKthRoot`
///  \sa `FieldWithRootOf`
///  \sa `AlgebraicStructureTraits`
class EuclideanRing {
public:

}; /* concept EuclideanRing */
/// @}
/// @} 

                   
  

