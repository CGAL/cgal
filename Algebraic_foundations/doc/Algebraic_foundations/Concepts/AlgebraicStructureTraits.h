/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup concepts Concepts
/// @{

 
///  
///  A model of `AlgebraicStructureTraits` reflects the algebraic structure
///  of an associated type `Type`.
///  Depending on the concepts that `Type` fulfills, 
///  it contains various functors and descriptive tags. 
///  Moreover it gives access to the several possible 
///  algebraic operations within that structure.
///  
///  A model of `AlgebraicStructureTraits` is supposed to provide:
///  
///  Functors 
///  -------------- 
///  
///  In case a functor is not provided, it is set to `CGAL::Null_functor`.
///  
///  \sa `IntegralDomainWithoutDivision`
///  \sa `IntegralDomain`
///  \sa `UniqueFactorizationDomain`
///  \sa `EuclideanRing`
///  \sa `Field`
///  \sa `FieldWithSqrt`
///  \sa `FieldWithKthRoot`
///  \sa `FieldWithRootOf`
///  \sa `CGAL::Integral_domain_without_division_tag`
///  \sa `CGAL::Integral_domain_tag`
///  \sa `CGAL::Unique_factorization_domain_tag`
///  \sa `CGAL::Euclidean_ring_tag`
///  \sa `CGAL::Field_tag`
///  \sa `CGAL::Field_with_sqrt_tag`
///  \sa `CGAL::Field_with_kth_root_tag`
///  \sa `CGAL::Field_with_root_of_tag`
///  \hasModel `CGAL::Algebraic_structure_traits`
class AlgebraicStructureTraits {
public:

/// \name Types
/// @{
/*!
 The associated type.
*/
typedef Hidden_type Type;
/// @}

/// \name Types
/// @{
/*!
  Tag indicating the algebraic structure of the associated type. 

  <TABLE CELLSPACING=5 >
  <TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR>
  <TR>
    <TD class="math" ALIGN=LEFT NOWRAP>
Tag is: 
    <TD class="math" ALIGN=LEFT NOWRAP>
`Type` is model of:
  <TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR>
  <TR>
    <TD class="math" ALIGN=LEFT NOWRAP>
`CGAL::Null_tag`                       
    <TD class="math" ALIGN=LEFT NOWRAP>
no algebraic concept
  <TR>
    <TD class="math" ALIGN=LEFT NOWRAP>
`CGAL::Integral_domain_without_division_tag`
    <TD class="math" ALIGN=LEFT NOWRAP>
`IntegralDomainWithoutDivision`
  <TR>
    <TD class="math" ALIGN=LEFT NOWRAP>
`CGAL::Integral_domain_tag`            
    <TD class="math" ALIGN=LEFT NOWRAP>
`IntegralDomain`
  <TR>
    <TD class="math" ALIGN=LEFT NOWRAP>
`CGAL::Unique_factorization_domain_tag`                   
    <TD class="math" ALIGN=LEFT NOWRAP>
`UniqueFactorizationDomain`
  <TR>
    <TD class="math" ALIGN=LEFT NOWRAP>
`CGAL::Euclidean_ring_tag`             
    <TD class="math" ALIGN=LEFT NOWRAP>
`EuclideanRing`
  <TR>
    <TD class="math" ALIGN=LEFT NOWRAP>
`CGAL::Field_tag`                      
    <TD class="math" ALIGN=LEFT NOWRAP>
`Field`
  <TR>
    <TD class="math" ALIGN=LEFT NOWRAP>
`CGAL::Field_with_sqrt_tag`            
    <TD class="math" ALIGN=LEFT NOWRAP>
`FieldWithSqrt`
  <TR>
    <TD class="math" ALIGN=LEFT NOWRAP>
`CGAL::Field_with_kth_root_tag`            
    <TD class="math" ALIGN=LEFT NOWRAP>
`FieldWithKthRoot`
  <TR>
    <TD class="math" ALIGN=LEFT NOWRAP>
`CGAL::Field_with_root_of_tag`            
    <TD class="math" ALIGN=LEFT NOWRAP>
`FieldWithRootOf`
  <TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR>
</TABLE>
        
        
*/
typedef Hidden_type Algebraic_category;
/// @}

/// \name Types
/// @{
/*!
  Tag indicating whether `Type` is exact. 
          This is either `CGAL::Tag_true` or `CGAL::Tag_false`.
          An algebraic structure is considered exact, if all operations 
          required by its concept are computed such that a comparison 
          of two algebraic expressions is always correct.
          The exactness covers only those operations that are required by 
          the algebraic structure concept. 
          e.g. an exact `Field` may have a `Sqrt` functor that 
          is not exact.    
        
*/
typedef Hidden_type Is_exact;
/// @}

/// \name Types
/// @{
/*!
  Tag indicating whether `Type` is numerical sensitive. 
          This is either `CGAL::Tag_true` or `CGAL::Tag_false`.
          An algebraic structure is  considered as numerically sensitive,
          if the performance of the type is sensitive to the condition 
          number of an algorithm.
        
*/
typedef Hidden_type Is_numerical_sensitive;
/// @}

/// \name Types
/// @{
/*!
  This type specifies the return type of the predicates provided
          by this traits. The type must be convertible to `bool` and 
          typically the type indeed maps to `bool`. However, there are also 
          cases such as interval arithmetic, in which it is `Uncertain<bool>` 
          or some similar type. 
        
*/
typedef Hidden_type Boolean;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::IsZero`.
Required by the concept `IntegralDomainWithoutDivision`.
In case `Type` is also model of `RealEmbeddable` this is a 
model of `RealEmbeddableTraits::IsZero`. 

*/
typedef Hidden_type Is_zero;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::IsOne`.
Required by the concept `IntegralDomainWithoutDivision`.

*/
typedef Hidden_type Is_one;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::Square`.
Required by the concept `IntegralDomainWithoutDivision`.

*/
typedef Hidden_type Square;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::Simplify`.
Required by the concept `IntegralDomainWithoutDivision`.

*/
typedef Hidden_type Simplify;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::UnitPart`.
Required by the concept `IntegralDomainWithoutDivision`.

*/
typedef Hidden_type Unit_part;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::IntegralDivision`.
Required by the concept `IntegralDomain`.

*/
typedef Hidden_type Integral_division;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::Divides`.
Required by the concept `IntegralDomain`.

*/
typedef Hidden_type Divides;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::IsSquare`.
Required by the concept `IntegralDomainWithoutDivision`.

*/
typedef Hidden_type Is_square;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::Gcd`.
Required by the concept `UniqueFactorizationDomain`.

*/
typedef Hidden_type Gcd;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::Mod`.
Required by the concept `EuclideanRing`.

*/
typedef Hidden_type Mod;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::Div`.
Required by the concept `EuclideanRing`.

*/
typedef Hidden_type Div;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::DivMod`.
Required by the concept `EuclideanRing`.

*/
typedef Hidden_type Div_mod;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::Inverse`.
Required by the concept `Field`.

*/
typedef Hidden_type Inverse;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::Sqrt`.
Required by the concept `FieldWithSqrt`.

*/
typedef Hidden_type Sqrt;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::KthRoot`.
Required by the concept `FieldWithKthRoot`.

*/
typedef Hidden_type Kth_root;
/// @}

/// \name Functors
/// @{
/*!
  
A model of `AlgebraicStructureTraits::RootOf`.
Required by the concept `FieldWithRootOf`.

*/
typedef Hidden_type Root_of;
/// @}

}; /* concept AlgebraicStructureTraits */
/// @}
/// @} 

                   
  

