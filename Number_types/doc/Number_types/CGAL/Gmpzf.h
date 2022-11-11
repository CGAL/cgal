
namespace CGAL {

/*!
\ingroup nt_gmp

An object of the class `Gmpzf`  is a multiple-precision floating-point number which can represent
numbers of the form \f$ m*2^e\f$, where \f$ m\f$ is an arbitrary precision integer
based on the \gmp library, and \f$ e\f$
is of type `long`. This type can be considered exact, even if the
exponent is not a multiple-precision number. This number type offers
functionality very similar to `MP_Float` but is generally faster.

\cgalModels `EuclideanRing`
\cgalModels `RealEmbeddable`

\cgalHeading{Implementation}

The significand \f$ m\f$ of a `Gmpzf` is a `Gmpz` and is reference
counted. The exponent \f$ e\f$ of a `Gmpzf` is a `long`.

*/

class Gmpzf {
public:

/// \name Creation
/// @{

/*!
creates a `Gmpzf` initialized with `0`.
*/
Gmpzf();

/*!
creates a `Gmpzf` initialized with `i`.
*/
Gmpzf(int i);

/*!
creates a `Gmpzf` initialized with `l`.
*/
Gmpzf(long int l);

/*!
creates a `Gmpzf` initialized with `i`.
*/
Gmpzf(const Gmpz& i);

/*!
creates a `Gmpzf` initialized with `f`.
*/
Gmpzf(const Gmpfr& f);

/*!
creates a `Gmpzf` initialized with `d`.
*/
Gmpzf(double d);

/// @}

}; /* end Gmpzf */

/*!
writes a double approximation of `f` to the ostream `out`.
\relates Gmpzf
*/
std::ostream& operator<<(std::ostream& out, const Gmpzf& f);

/*!
writes an exact representation of `f` to the ostream `out`.
\relates Gmpzf
*/
std::ostream& print (std::ostream& out, const Gmpzf& f);

/*!
reads a `double` from `in`, then converts it to a `Gmpzf`.
\relates Gmpzf
*/
std::istream& operator>>(std::istream& in, Gmpzf& f);

} /* end namespace CGAL */

