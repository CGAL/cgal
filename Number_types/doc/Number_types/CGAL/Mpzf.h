
namespace CGAL {

/*!
\ingroup nt_gmp

An object of the class `Mpzf` is a multiple-precision floating-point
number which can represent numbers of the form \f$ m*2^e\f$, where \f$
m\f$ is an arbitrary precision integer based on the \gmp
library, and \f$ e\f$ is of type `int`. This
type can be considered exact, even if the exponent is not a
multiple-precision number. A `Mpzf` constructed from an integer or a
normalized finite floating point number represents exactly that number.
This number type offers functionality very
similar to `MP_Float` and `Gmpzf` but is faster.

\cgalModels `IntegralDomainWithoutDivision`
\cgalModels `RealEmbeddable`

\cgalHeading{Implementation}

This class is only available on platforms on which \gmp uses 64 bit limbs and the endianness is
either big-endian or little-endian. The macro `CGAL_HAS_MPZF` will be defined
when including `CGAL/Mpzf.h` if this is true.  This class makes the assumption that the
representation of a `double` in memory follows IEEE 754.

Currently, an `Mpzf` contains an array of a few limbs, in which it
stores the data as long as it fits. If it does not fit, it dynamically
allocates memory with new, and de-allocates it when it is done. Code to
recycle the allocated memory (per thread) is included but disabled at
the moment.
*/

struct Mpzf {

/// \name Creation
/// @{

/*!
creates a `Mpzf` initialized with `0`.
*/
Mpzf();

/*!
creates a `Mpzf` initialized with `i`.
*/
Mpzf(int i);

/*!
creates a `Mpzf` initialized with `l`.
*/
Mpzf(long int l);

/*!
creates a `Mpzf` initialized with `i`.
*/
Mpzf(const Gmpz& i);

/*!
creates a `Mpzf` initialized with `i`.
*/
Mpzf(const mpz_class& i);

/*!
creates a `Mpzf` initialized with `d`.
*/
Mpzf(double d);

/// @}

/// \name Conversion
/// @{

/*!
creates a `Gmpq` initialized with `*this`.
*/
explicit operator Gmpq() const;

/*!
creates a `mpq_class` initialized with `*this`.
*/
explicit operator mpq_class() const;

/// @}

}; /* end Mpzf */

/*!
writes a double approximation of `f` to the ostream `out`.
\relates Mpzf
*/
std::ostream& operator<<(std::ostream& out, const Mpzf& f);

/*!
reads a `double` from `in`, then converts it to a `Mpzf`.
\relates Mpzf
*/
std::istream& operator>>(std::istream& in, Mpzf& f);

} /* end namespace CGAL */

