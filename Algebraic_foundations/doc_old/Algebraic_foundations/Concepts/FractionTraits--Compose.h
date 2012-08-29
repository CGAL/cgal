/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  `AdaptableBinaryFunction`, returns the fraction of its arguments.
///  \refines ::AdaptableBinaryFunction
///  \sa `Fraction`
///  \sa `FractionTraits`
///  \sa `FractionTraits::Decompose`
///  \sa `FractionTraits::CommonFactor`
class FractionTraits::Compose {
public:

/// \name Types
/// @{
/*!
 
*/
typedef FractionTraits::Type    result_type;
/// @}

/// \name Types
/// @{
/*!
 
*/
typedef FractionTraits::Numerator_type   first_argument_type;
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
 return the fraction \f$n/d\f$. 
*/
result_type operator()(first_argument_type n, second_argument_type d);
/// @}

}; /* concept FractionTraits::Compose */
/// @}
/// @} 

                   
  

