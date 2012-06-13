///
namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{
///
 
///  
///  The function `simplest_rational_in_interval` computes the simplest rational number in an
///  interval of two `double` values.
///  
///  Implementation 
///  -------------- 
///  
///  See Knuth, "Seminumerical algorithms", page 654, answer to exercise
///  4.53-39.
/// computes the rational number with the smallest denominator in the  interval `[d1,d2].` 
Rational simplest_rational_in_interval(double d1, double d2);
/// @}
} 

