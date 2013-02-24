
namespace CGAL {

/*!
\ingroup nt_gmp

An object of the class `mpzf` is a multiple-precision floating-point
number which can represent numbers of the form \f$ m*2^e\f$, where \f$
m\f$ is an arbitrary precision integer based on the <span
class="textsc">Gmp</span> library, and \f$ e\f$ is of type `int`. This
type can be considered exact, even if the exponent is not a
multiple-precision number. This number type offers functionality very
similar to `MP_Float` and `Gmpzf` but is faster.

\cgalModels `IntegralDomainWithoutDivision`
\cgalModels `RealEmbeddable`

\cgalHeading{Implementation}

This class is only available on platforms on which <span
class="textsc">Gmp</span> uses 64 bits limbs and the endianness is
either big-endian or little-endian.  If this is true, the macro
`CGAL_HAS_MPZF` is defined.  This class makes the assumption that the
representation of a `double` in memory follows IEEE 754.

Currently, an `mpzf` contains an array of a few limbs, in which it
stores the data as long as it fits. If it does not fit, it dynamically
allocates memory with new, and de-allocates it when it is done. Code to
recycle the allocated memory (per thread) is included but disabled at
the moment.
*/

class mpzf {
public:

/// \name Creation
/// @{

/*!
creates a `mpzf` initialized with `0`.
*/
mpzf();

/*!
creates a `mpzf` initialized with `i`.
*/
mpzf(int i);

/*!
creates a `mpzf` initialized with `l`.
*/
mpzf(long int l);

/*!
creates a `mpzf` initialized with `i`.
*/
mpzf(const Gmpz& i);

/*!
creates a `mpzf` initialized with `i`.
*/
mpzf(const mpz_class& i);

/*!
creates a `mpzf` initialized with `d`.
*/
mpzf(double d);

/// @}

/// \name Conversion
/// @{

/*!
creates a `Gmpq` initialized with `*this`.
*/
operator Gmpq() const;

/*!
creates a `mpq_class` initialized with `*this`.
*/
operator mpq_class() const;

/// @}

}; /* end mpzf */

/*!
writes a double approximation of `f` to the ostream `out`.
\relates mpzf
*/
std::ostream& operator<<(std::ostream& out, const mpzf& f);

/*!
reads a `double` from `in`, then converts it to a `mpzf`.
\relates mpzf
*/
std::istream& operator>>(std::istream& in, mpzf& f);

} /* end namespace CGAL */

