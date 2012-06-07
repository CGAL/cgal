/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup concepts Concepts
/// @{

 
///  
///  This `AdaptableUnaryFunction` computes the sign of a real embeddable number.
///  \refines `AdaptableUnaryFunction`
///  \sa `RealEmbeddableTraits`
class RealEmbeddableTraits::Sgn {
public:

/// \name Types
/// @{
/*!
  Type convertible to `CGAL::Sign`.
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
 Computes the sign of \f$x\f$. 
*/
result_type operator()(argument_type x);
/// @}

}; /* concept RealEmbeddableTraits::Sgn */
/// @}
/// @} 

                   
  

