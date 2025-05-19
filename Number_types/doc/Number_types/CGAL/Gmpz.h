
namespace CGAL {

/*!
\ingroup nt_gmp

An object of the class `Gmpz` is an arbitrary precision integer
based on the \gmp Library.

\cgalModels{EuclideanRing,RealEmbeddable}

\cgalHeading{Implementation}

`Gmpz`s are reference counted.

*/

class Gmpz {
public:

/// \name Creation
/// @{

/*!
creates an uninitialized multiple precision integer `z`.
*/
Gmpz();

/*!
creates a multiple-precision integer initialized with `i`.
*/
Gmpz(int i);

/*!
creates a multiple-precision integer initialized with
the integral part of `d`.
*/
Gmpz(double d);

/// @}

/// \name Operations
/// @{

/*!
prefix increment.
*/
Gmpz & operator++();

/*!
postfix increment.
*/
Gmpz operator++(int);

/*!
prefix decrement.
*/
Gmpz & operator--();

/*!
postfix decrement.
*/
Gmpz operator--(int);

/*!
rightshift by `i`, where `i >= 0`.
*/
Gmpz & operator>>=(const long& i);

/*!
leftshift by `i`, where  `i >= 0`.
*/
Gmpz & operator<<=(const long& i);

/*!
bitwise AND.
*/
Gmpz & operator&=(const Gmpz& b);

/*!
bitwise IOR.
*/
Gmpz & operator|=(const Gmpz& b);

/*!
bitwise XOR.
*/
Gmpz & operator^=(const Gmpz& b);

/*!
Returns the sign of `z`.
*/
Sign sign() const;

/*!
Returns the bit-size (that is, the number of bits needed to
represent the mantissa) of `z`.
*/
size_t bit_size() const;

/*!
Returns the size in limbs of `z`. A limb is the type used by
\gmp to represent the integer (usually `long`).
*/
size_t size() const;

/*!
Returns the approximate number of decimal digits needed to
represent `z`. Approximate means either a correct result, either
the correct result plus one.
*/
size_t approximate_decimal_length() const;

/*!
Returns a double approximation of `z`. The integer is truncated
if needed. If the exponent of the conversion is too big, the result
is system dependent (returning infinity where it is supported).
*/
double to_double() const;

/// @}

}; /* end Gmpz */

/*!
rightshift by `i`.
\relates Gmpz
*/
Gmpz operator>>(const Gmpz& a, unsigned long i);

/*!
leftshift by `i`.
\relates Gmpz
*/
Gmpz operator<<(const Gmpz& a, unsigned long i);


/*!
bitwise AND.
\relates Gmpz
*/
Gmpz operator&(const Gmpz& a, const Gmpz& b);

/*!
bitwise IOR.
\relates Gmpz
*/
Gmpz operator|(const Gmpz& a, const Gmpz& b);

/*!
bitwise XOR.
\relates Gmpz
*/
Gmpz operator^(const Gmpz& a, const Gmpz& b);

/*!
writes `z` to the ostream `out`.
\relates Gmpz
*/
std::ostream& operator<<(std::ostream& out, const Gmpz& z);

/*!
reads an integer from `in`, then converts it to a `Gmpz`.
\relates Gmpz
*/
std::istream& operator>>(std::istream& in, Gmpz& z);


} /* end namespace CGAL */
