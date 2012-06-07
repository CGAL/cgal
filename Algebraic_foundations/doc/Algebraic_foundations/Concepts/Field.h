/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup concepts Concepts
/// @{

 
///  
///  A model of `Field` is an `IntegralDomain` in which every non-zero element
///  has a multiplicative inverse. 
///  Thus, one can divide by any non-zero element. 
///  Hence division is defined for any divisor != 0. 
///  For a Field, we require this division operation to be available through 
///  operators / and /=.
///
///  Moreover, `CGAL::Algebraic_structure_traits` is a model of 
///  `AlgebraicStructureTraits` providing:
///  - `CGAL::Algebraic_structure_traits::Algebraic_type` derived from `Field_tag` 
///  - `CGAL::Algebraic_structure_traits
///
///  \refines `IntegralDomain`
///  \sa `IntegralDomainWithoutDivision`
///  \sa `IntegralDomain`
///  \sa `UniqueFactorizationDomain`
///  \sa `EuclideanRing`
///  \sa `Field`
///  \sa `FieldWithSqrt`
///  \sa `FieldWithKthRoot`
///  \sa `FieldWithRootOf`
///  \sa `AlgebraicStructureTraits`
class Field {
public:

/*!
 
*/
Field operator/=(const Field &b);

}; /* concept Field */

/// @}
/// @} 

/// \relates Field
/// 
Field operator/(const Field &a, const Field &b);

                   
  

