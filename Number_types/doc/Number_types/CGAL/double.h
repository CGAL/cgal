/// \addtogroup PkgNumberTypes Number Types
/// @{

 
///  The fundamental type `double` is an `RealEmbeddable` 
///  `Field`. Due to rounding errors and overflow `double` is 
///  considered as not exact.
///  
///  \models ::FieldWithSqrt
///  \models ::RealEmbeddable
class double {
public:

 
}; /* class double */

namespace CGAL {
/// \relates double
/// Determines whether the argument represents a value in \f$R\f$. 
bool is_finite(double x);

} // CGAL

/// @}
