namespace CGAL {
/// \addtogroup PkgNumberTypes Number Types
/// @{

 
///  
///  An object of the class `Gmpz` is an arbitrary precision integer
///  based on the <span class="textsc">Gnu</span> Multiple Precision Arithmetic Library.
///  
///  \models ::EuclideanRing
///  \models ::RealEmbeddable
///
///  Implementation 
///  -------------- 
///  `Gmpz`s are reference counted.
///   
///  }; /* class Gmpz */
class Gmpz {
public:
/// \name Creation
/// @{
/*!
 creates an uninitialized multiple precision integer `z`.
*/
Gmpz();
/// @}

/// \name Creation
/// @{
/*!
 creates a multiple-precision integer initialized with              `i`.
*/
  Gmpz(int i);
/// @}

/// \name Creation
/// @{
/*!
 creates a multiple-precision integer initialized with              the integral part of `d`.
*/
  Gmpz(double d);
/// @}

/// \name Operations
/// @{
/*!
 prefix increment.
*/
Gmpz & operator++();
/// @}

/// \name Operations
/// @{
/*!
 postfix increment.
*/
Gmpz   operator++(int);
/// @}

/// \name Operations
/// @{
/*!
 prefix decrement.
*/
Gmpz & operator--();
/// @}

/// \name Operations
/// @{
/*!
 postfix decrement.
*/
Gmpz   operator--(int);
/// @}

/// \name Operations
/// @{
/*!
 rightshift by i, where \f$i>=0\f$.
*/
Gmpz & operator>>=(const long& i);
/// @}

/// \name Operations
/// @{
/*!
 leftshift by i, where \f$i>=0\f$.
*/
Gmpz & operator<<=(const long& i);
/// @}

/// \name Operations
/// @{
/*!
 bitwise AND.
*/
Gmpz & operator&=(const Gmpz& b);
/// @}

/// \name Operations
/// @{
/*!
 bitwise IOR.
*/
Gmpz & operator|=(const Gmpz& b);
/// @}

/// \name Operations
/// @{
/*!
 bitwise XOR.
*/
Gmpz & operator^=(const Gmpz& b);
/// @}

/// \name Operations
/// @{
/*!
 Returns the sign of `z`.
*/
Sign sign() const;
/// @}

/// \name Operations
/// @{
/*!
 Returns the bit-size (that is, the number of bits needed to         represent the mantissa) of `z`.
*/
size_t bit_size() const;
/// @}

/// \name Operations
/// @{
/*!
 Returns the size in limbs of `z`. A limb is the type used by         <span class="textsc">Gmp</span> to represent the integer (usually `long`).
*/
size_t size() const;
/// @}

/// \name Operations
/// @{
/*!
 Returns the approximate number of decimal digits needed to         represent `z`. Approximate means either a correct result, either         the correct result plus one.
*/
size_t approximate_decimal_length() const;
/// @}

/// \name Operations
/// @{
/*!
 Returns a double approximation of `z`. The integer is truncated         if needed. If the exponent of the conversion is too big, the result         is system dependent (returning infinity where it is supported).
*/
double to_double() const;
/// @}

}; /* class Gmpz */

///
/// \relates Gmpz
/// rightshift by \f$i\f$.
Gmpz operator>>(const Gmpz& a, unsigned long i);
///
/// \relates Gmpz
/// leftshift by \f$i\f$.
Gmpz operator<<(const Gmpz& a, unsigned long i);

///
/// \relates Gmpz
/// writes `z` to the ostream `out`.
std::ostream& operator<<(std::ostream& out, const Gmpz& z);

///
/// \relates Gmpz
/// reads an integer from `in`, then converts it to a `Gmpz`.
std::istream& operator>>(std::istream& in, Gmpz& z);

///
/// \relates Gmpz
/// bitwise AND.
Gmpz operator&(const Gmpz& a, const Gmpz& b);
///
/// \relates Gmpz
/// bitwise IOR.
Gmpz operator|(const Gmpz& a, const Gmpz& b);
///
/// \relates Gmpz
/// bitwise XOR.
Gmpz operator^(const Gmpz& a, const Gmpz& b);

/// @}
} // namespace CGAL

                   
  

