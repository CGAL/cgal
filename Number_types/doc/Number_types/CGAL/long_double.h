/// \addtogroup PkgNumberTypes Number Types
/// @{
 
///  The fundamental type `long double` is an `RealEmbeddable` 
///  `FieldWithSqrt`. Due to rounding errors and overflow `long double` is 
///  considered as not exact.
///  
///  \models ::FieldWithSqrt
///  \models ::RealEmbeddable
class long_double {
public:

}; /* class long double */

/// @}

namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{

///
/// \relates long_double
/// Determines whether the argument represents a value in \f$R\f$. 
bool is_finite(long double x);
/// @}

} // namespace CGAL

                   
  

