/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup concepts Concepts
/// @{

 
///  
///  `AdaptableBinaryFunction`, finds great common factor of denominators.
///  This can be considered as a relaxed version of `AlgebraicStructureTraits::Gcd`, 
///  this is needed because it is not guaranteed that `FractionTraits::Denominator_type` is a model of 
///  `UniqueFactorizationDomain`.
///  \refines `AdaptableBinaryFunction`
///  \sa `Fraction`
///  \sa `FractionTraits`
///  \sa `FractionTraits::Decompose`
///  \sa `FractionTraits::Compose`
///  \sa `AlgebraicStructureTraits::Gcd`
class FractionTraits::CommonFactor {
public:

/// \name Types
/// @{
/*!
 
*/
typedef FractionTraits::Denominator_type result_type;
/// @}

/// \name Types
/// @{
/*!
 
*/
typedef FractionTraits::Denominator_type first_argument_type;
/// @}

/// \name Types
/// @{
/*!
 
*/
typedef FractionTraits::Denominator_type second_argument_type;
/// @}

/// \name Operations
/// @{
/*!
 return a great common factor of \f$d1\f$ and \f$d2\f$.          Note: <TT>operator()(0,0) = 0</TT> 
*/
result_type operator()(first_argument_type d1, second_argument_type d2);
/// @}

}; /* concept FractionTraits::CommonFactor */
/// @}
/// @} 

                   
  

