
namespace CGAL {

/*!
\ingroup nt_gmp

An object of the class `Gmpq` is an arbitrary precision rational
number based on the \gmp library.

\cgalModels `Field`
\cgalModels `RealEmbeddable`
\cgalModels `Fraction`

\cgalHeading{Implementation}

`Gmpq`s are reference counted.

*/

class Gmpq {
public:

/// \name Creation
/// @{

/*!
creates an uninitialized `Gmpq` `q`.
*/
Gmpq();

/*!
creates a `Gmpq` initialized with `i`.
*/
Gmpq(int i);

/*!
creates a `Gmpq` initialized with `n`.
*/
Gmpq(Gmpz n);

/*!
creates a `Gmpq` initialized with `f`.
*/
Gmpq(Gmpfr f);

/*!
creates a `Gmpq` initialized with `n/d`.
*/
Gmpq(int n, int d);

/*!
creates a `Gmpq` initialized with `n/d`.
*/
Gmpq(signed long n, unsigned long d);

/*!
creates a `Gmpq` initialized with `n/d`.
*/
Gmpq(unsigned long n, unsigned long d);

/*!
creates a `Gmpq` initialized with `n/d`.
*/
Gmpq(Gmpz n, Gmpz d);

/*!
creates a `Gmpq` initialized with `d`.
*/
Gmpq(double d);

/*!
creates a `Gmpq` initialized with `str`, which can
be an integer like "41" or a fraction like "41/152". White
space is allowed in the string, and ignored.
*/
Gmpq(const std::string& str);

/*!
creates a `Gmpq` initialized with `str` in base
`base`, which is an integer between 2 and 62. White space
in the string is ignored.
*/
Gmpq(const std::string& str, int base);

/// @}

/// \name Operations
/// There are two access functions, namely to the numerator and the
/// denominator of a rational. Note that these values are not uniquely
/// defined. It is guaranteed that `q.numerator()` and
/// `q.denominator()` return values `nt_num` and `nt_den` such that `q
/// = nt_num/nt_den`, only if `q.numerator()` and `q.denominator()`
/// are called consecutively wrt. `q`, i.e.\ `q` is not involved in any
/// other operation between these calls.
/// @{

/*!
returns the numerator of `q`.
*/
Gmpz numerator() const;

/*!
returns the denominator of `q`.
*/
Gmpz denominator() const;

/// @}

}; /* end Gmpq */

/*!
writes `q` to the ostream `out`, in the form `n/d`.
\relates Gmpq
*/
std::ostream& operator<<(std::ostream& out, const Gmpq& q);

/*!
reads a number from `in`, then converts it to a
`Gmpq`. The number may be an integer, a rational number in
the form `n/d`, or a floating-point number.
\relates Gmpq
*/
std::istream& operator>>(std::istream& in, Gmpq& q);

} /* end namespace CGAL */
