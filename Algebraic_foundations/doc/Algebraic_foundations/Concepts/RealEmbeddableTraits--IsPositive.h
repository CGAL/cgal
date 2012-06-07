/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup concepts Concepts
/// @{

 
///  
///  `AdaptableUnaryFunction`, returns true in case the argument is positive.
///  \refines `AdaptableUnaryFunction`
///  \sa `RealEmbeddableTraits`
class RealEmbeddableTraits::IsPositive {
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
 returns true in case \f$x\f$ is positive. 
*/
result_type operator()(argument_type x);
/// @}

}; /* concept RealEmbeddableTraits::IsPositive */
/// @}
/// @} 

                   
  

