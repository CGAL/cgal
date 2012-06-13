namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{

 
///  An object of the class `MP_Float` is able to represent a floating point
///  value with arbitrary precision.  This number type has the property that
///  additions, subtractions and multiplications are computed exactly, as well as
///  the construction from `float`, `double` and `long double`.
///  Division and square root are not enabled by default since \cgal release 3.2,
///  since they are computed approximately.  We suggest that you use
///  rationals like `Quotient<MP_Float>` when you need exact divisions.
///  Note on the implementation : although the mantissa length is basically only
///  limited by the available memory, the exponent is currently represented by a
///  (integral valued) `double`, which can overflow in some circumstances.  We
///  plan to also have a multiprecision exponent to fix this issue.
///  
///  \models ::EuclideanRing
///  \models ::RealEmbeddable
///
///  Implementation 
///  -------------- 
///  The implementation of `MP_Float` is simple but provides a quadratic
///  complexity for multiplications.  This can be a problem for large operands.
///  For faster implementations of the same functionality with large integral
///  values, you may want to consider using `GMP` or `LEDA` instead.
///   
class MP_Float {
public:

/// \name Creation
/// @{
/*!
 introduces an uninitialized variable `m`.
*/
MP_Float();
/// @}

/// \name Creation
/// @{
/*!
 copy constructor.
*/
MP_Float(const MP_Float &);
/// @}

/// \name Creation
/// @{
/*!
 introduces the integral value i.
*/
  MP_Float(int i);
/// @}

/// \name Creation
/// @{
/*!
 introduces the floating point value d (exact conversion).
*/
  MP_Float(float d);
/// @}

/// \name Creation
/// @{
/*!
 introduces the floating point value d (exact conversion).
*/
  MP_Float(double d);
/// @}

/// \name Creation
/// @{
/*!
 introduces the floating point value d (exact conversion).
*/
  MP_Float(long double d);
/// @}

}; /* class MP_Float */

/// \relates MP_Float
/// writes a double approximation of `m` to the ostream `out`.
std::ostream& operator<<(std::ostream& out, const MP_Float& m);

/// \relates MP_Float
/// reads a `double` from `in`, then converts it to an `MP_Float`.
std::istream& operator>>(std::istream& in, MP_Float& m);

/// \relates MP_Float
/// computes an approximation of the division by converting the operands to `double`, performing the division on `double`, and converting back to `MP_Float`.
MP_Float approximate_division(const MP_Float &a, const MP_Float &b);

/// \relates MP_Float
/// computes an approximation of the square root by converting the operand to `double`, performing the square root on `double`, and converting back to `MP_Float`.
MP_Float approximate_sqrt(const MP_Float &a);

/// @}
} // namespace CGAL

                   
  

