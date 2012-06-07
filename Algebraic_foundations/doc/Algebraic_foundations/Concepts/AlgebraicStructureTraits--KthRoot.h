/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup concepts Concepts
/// @{

 
///  
///  `AdaptableBinaryFunction` providing the k-th root.
///  \refines `AdaptableBinaryFunction`
///  \sa `FieldWithRootOf`
///  \sa `AlgebraicStructureTraits`
class AlgebraicStructureTraits::KthRoot {
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
  Is int.
*/
typedef Hidden_type first_argument;
/// @}

/// \name Types
/// @{
/*!
  Is `AlgebraicStructureTraits::Type`.
*/
typedef Hidden_type second_argument;
/// @}

/// \name Operations
/// @{
/*!
 returns the \f$k\f$-th root of \f$x\f$.          \pre  \f$k   \geq 1\f$           
*/
result_type operator()(int k, second_argument_type x);
/// @}

}; /* concept AlgebraicStructureTraits::KthRoot */
/// @}
/// @} 

                   
  

