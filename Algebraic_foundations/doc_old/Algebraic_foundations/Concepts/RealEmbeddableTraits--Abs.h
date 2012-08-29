/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  `AdaptableUnaryFunction` computes the absolute value of a number.
///  \refines ::AdaptableUnaryFunction
///  \sa `RealEmbeddableTraits`
class RealEmbeddableTraits::Abs {
public:

/// \name Types
/// @{
/*!
  Is `RealEmbeddableTraits::Type`.
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
 computes the absolute value of \f$x\f$. 
*/
result_type operator()(argument_type  x);
/// @}

}; /* concept RealEmbeddableTraits::Abs */
/// @}
/// @} 

                   
  

