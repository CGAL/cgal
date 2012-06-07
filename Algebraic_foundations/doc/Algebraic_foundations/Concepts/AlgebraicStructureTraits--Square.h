/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup concepts Concepts
/// @{

 
///  
///  `AdaptableUnaryFunction`, computing the square of the argument.
///  \refines `AdaptableUnaryFunction`
///  \sa `AlgebraicStructureTraits`
class AlgebraicStructureTraits::Square {
public:

/// \name Types
/// @{
/*!
  Is `AlgebraicStructureTraits::Type`.
*/
typedef Hidden_type result_type;
/// @}

/// \name Types
/// @{
/*!
  Is `AlgebraicStructureTraits::Type`.
*/
typedef Hidden_type argument_type;
/// @}

/// \name Operations
/// @{
/*!
 returns the square of \f$x\f$.
*/
result_type operator()(argument_type x);
/// @}

}; /* concept AlgebraicStructureTraits::Square */
/// @}
/// @} 

                   
  

