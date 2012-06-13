/// \addtogroup PkgNumberTypes Number Types
/// @{

 
///  The fundamental type `float` is an `RealEmbeddable` 
///  `FieldWithSqrt`. Due to rounding errors and overflow `float` is 
///  considered as not exact.
///  
///  \models ::FieldWithSqrt
///  \models ::RealEmbeddable
class float {
public:
 
}; /* class float */

/// @}

namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{

/// \relates float
/// Determines whether the argument represents a value in \f$R\f$. 
bool is_finite(float x);

/// @}
}
