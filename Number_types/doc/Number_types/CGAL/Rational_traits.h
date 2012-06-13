namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{

 
///  
///  The class `Rational_traits` can be used to determine the type of the numerator
///  and denominator of a rational number type as `Quotient`, `Gmpq`,
///   `mpq_class` or `leda_rational`.
///  
///  
///  
template< class NT >
class Rational_traits {
public:

/// \name Types
/// @{
/*!
 the type of the numerator and denominator.
*/
typedef Hidden_type RT;
/// @}

/// \name Operations
/// @{
/*!
 returns the numerator of `r`.
*/
RT numerator   (const NT & r) const;
/// @}

/// \name Operations
/// @{
/*!
 returns the denominator of `r`.
*/
RT denominator (const NT & r) const;
/// @}

/// \name Operations
/// @{
/*!
 constructs a rational number.
*/
NT make_rational(const RT & n, const RT & d) const;
/// @}

/// \name Operations
/// @{
/*!
 constructs a rational number.
*/
NT make_rational(const NT & n, const NT & d) const;
/// @}

 
}; /* class Rational_traits */
/// @}
} // namespace CGAL

                   
  

