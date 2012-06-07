/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup concepts Concepts
/// @{

 
///  
///  This is the most basic concept for algebraic structures considered within CGAL.
///  A model `IntegralDomainWithoutDivision`  represents an integral domain,
///  i.e. commutative ring with 0, 1, +, * and unity free of zero divisors. 
///  <B>Note:</B> A model is not required to offer the always well defined integral division.
///
///  It refines `Assignable`, `CopyConstructible`, `DefaultConstructible`
///  and `FromIntConstructible`.  
///  It refines `EqualityComparable`, where equality is defined w.r.t. 
///  the ring element being represented. 
///  The operators unary and binary plus +, unary and binary minus -, 
///  multiplication * and their compound forms +=, -=, *= are required and 
///  implement the respective ring operations.
///
///  Moreover, `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >` is a model of 
///  `AlgebraicStructureTraits` providing: 
///  - `CGAL::Algebraic_structure_traits::Algebraic_type` derived from `Integral_domain_without_division_tag`
///  - `CGAL::Algebraic_structure_traits::Is_zero`  
///  - `CGAL::Algebraic_structure_traits::Is_one`  
///  - `CGAL::Algebraic_structure_traits::Square`  
///  - `CGAL::Algebraic_structure_traits::Simplify` 
///  - `CGAL::Algebraic_structure_traits::Unit_part` 
///  
///  \refines `Assignable`
///  \refines `CopyConstructible`
///  \refines `DefaultConstructible`
///  \refines `EqualityComparable`
///  \refines `FromIntConstructible`
///
///  \sa `IntegralDomainWithoutDivision`
///  \sa `IntegralDomain`
///  \sa `UniqueFactorizationDomain`
///  \sa `EuclideanRing`
///  \sa `Field`
///  \sa `FieldWithSqrt`
///  \sa `FieldWithKthRoot`
///  \sa `FieldWithRootOf`
///  \sa `AlgebraicStructureTraits`
class IntegralDomainWithoutDivision {
public:
/// \name Operations
/// @{
/*!
 
*/
IntegralDomainWithoutDivision 
            operator+=(const IntegralDomainWithoutDivision &b);
/// @}

/// \name Operations
/// @{
/*!
 
*/
IntegralDomainWithoutDivision 
            operator-=(const IntegralDomainWithoutDivision &b);
/// @}

/// \name Operations
/// @{
/*!
 
*/
IntegralDomainWithoutDivision 
            operator*=(const IntegralDomainWithoutDivision &b);
/// @}

}; /* concept IntegralDomainWithoutDivision */

/// @}
/// @} 

///
/// \relates IntegralDomainWithoutDivision
/// unary plus 
IntegralDomainWithoutDivision 
            operator+(const IntegralDomainWithoutDivision &a);

///
/// \relates IntegralDomainWithoutDivision
/// unary minus
IntegralDomainWithoutDivision 
            operator-(const IntegralDomainWithoutDivision &a);

///
/// \relates IntegralDomainWithoutDivision
/// 
IntegralDomainWithoutDivision 
            operator+(const IntegralDomainWithoutDivision &a, 
                      const IntegralDomainWithoutDivision &b);

///
/// \relates IntegralDomainWithoutDivision
/// 
IntegralDomainWithoutDivision 
            operator-(const IntegralDomainWithoutDivision &a, 
                      const IntegralDomainWithoutDivision &b);

///
/// \relates IntegralDomainWithoutDivision
/// 
IntegralDomainWithoutDivision 
            operator*(const IntegralDomainWithoutDivision &a, 
                      const IntegralDomainWithoutDivision &b);


/// \relates IntegralDomainWithoutDivision
/// The `result_type` is convertible to `bool`. 
 result_type 
            operator==(const IntegralDomainWithoutDivision &a, 
                      const IntegralDomainWithoutDivision &b);

///
/// \relates IntegralDomainWithoutDivision
/// The `result_type` is convertible to `bool`. 
 result_type 
            operator!=(const IntegralDomainWithoutDivision &a, 
                      const IntegralDomainWithoutDivision &b);

                   
  

