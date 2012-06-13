/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  `AdaptableBinaryFunction` computes the integral quotient of division 
///  with remainder.
///  \refines ::AdaptableBinaryFunction
///  \sa `AlgebraicStructureTraits`
///  \sa `AlgebraicStructureTraits::Mod`
///  \sa `AlgebraicStructureTraits::DivMod`
class AlgebraicStructureTraits::Div {
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
 
*/
result_type operator()(first_argument_type  x, 
                                 second_argument_type y);
/// @}

/// \name Operations
/// @{
/*!
 This operator is defined if `NT1` and `NT2` are `ExplicitInteroperable`            with coercion type `AlgebraicStructureTraits::Type`. 
*/
template <class NT1, class NT2> result_type operator()(NT1  x, NT2  y);
/// @}

}; /* concept AlgebraicStructureTraits::Div */
/// @}
/// @} 

                   
  

