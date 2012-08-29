/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  `AdaptableUnaryFunction` providing the inverse element with 
///  respect to multiplication of a `Field`.
///  \refines ::AdaptableUnaryFunction
///  \sa `AlgebraicStructureTraits`
class AlgebraicStructureTraits::Inverse {
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
 returns  the inverse element of \f$x\f$ with respect to multiplication.            \pre \f$x  \neq 0\f$         
*/
result_type operator()(argument_type  x) const;
/// @}

}; /* concept AlgebraicStructureTraits::Inverse */
/// @}
/// @} 

                   
  

