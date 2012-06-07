/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup concepts Concepts
/// @{

 
///  
///  Functor decomposing a `Fraction` into its numerator and denominator.
///  
///  
///  
///  
///  \sa `Fraction`
///  \sa `FractionTraits`
///  \sa `FractionTraits::Compose`
///  \sa `FractionTraits::CommonFactor`
class FractionTraits::Decompose {
public:

/// \name Operations
/// @{
/*!
 decompose \f$f\f$ into numerator \f$n\f$ and denominator \f$d\f$. 
*/
void operator()( FractionTraits::Type               f,
                           FractionTraits::Numerator_type   & n,
                           FractionTraits::Denominator_type & d);
/// @}

}; /* concept FractionTraits::Decompose */
/// @}
/// @} 

                   
  

