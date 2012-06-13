/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  `AdaptableUnaryFunction`, returns true in case the argument is negative.
///  \refines ::AdaptableUnaryFunction
///  \sa `RealEmbeddableTraits`
class RealEmbeddableTraits::IsNegative {
public:

/// \name Types
/// @{
/*!
  Type convertible to `bool`.
*/
typedef Hidden_type result_type;
/// @}

/// \name Types
/// @{
/*!
  Is `RealEmbeddableTraits::Type`.
*/
typedef Hidden_type argument_type;
/// @}

/// \name Operations
/// @{
/*!
 returns true in case \f$x\f$ is negative. 
*/
result_type operator()(argument_type x);
/// @}

}; /* concept RealEmbeddableTraits::IsNegative */
/// @}
/// @} 

                   
  

