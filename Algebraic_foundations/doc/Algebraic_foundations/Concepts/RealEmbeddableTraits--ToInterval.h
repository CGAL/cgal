/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  `AdaptableUnaryFunction` computes for a given real embeddable 
///  number \f$x\f$ a double interval containing \f$x\f$. 
///  This interval is represented by `std::pair<double,double>`. 
///  
///  \refines ::AdaptableUnaryFunction
///  \sa `RealEmbeddableTraits`
class RealEmbeddableTraits::ToInterval {
public:

/// \name Types
/// @{
/*!
  Is `std::pair<double,double>`.
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
 computes a double interval containing \f$x\f$.
*/
result_type operator()(argument_type  x);
/// @}

}; /* concept RealEmbeddableTraits::ToInterval */
/// @}
/// @} 

                   
  

