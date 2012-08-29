/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  `AdaptableUnaryFunction`, returns true in case the argument is 0.
///  \refines ::AdaptableUnaryFunction
///  \sa `RealEmbeddableTraits`
///  \sa `AlgebraicStructureTraits::IsZero`
class RealEmbeddableTraits::IsZero {
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
 returns true in case \f$x\f$ is the zero element of the ring. 
*/
result_type operator()(argument_type x);
/// @}

}; /* concept RealEmbeddableTraits::IsZero */
/// @}
/// @} 

                   
  

