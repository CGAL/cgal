
namespace CGAL {

/*!
\ingroup nt_gmp

An object of the class `Gmpfr` is a fixed precision floating-point 
number, based on the <span class="textsc">Mpfr</span> library. This type is inexact, due to the fact 
that the mantissa of each number is represented by a fixed amount of bits 
(this amount is called <I>precision</I>). If an operation needs more bits 
than the precision of the result number, the results are rounded following 
different possible criteria (called <I>rounding modes</I>). 

Currently, the `Gmpfr` interface supports four rounding modes: round to nearest,
round toward zero, round down (or toward \f$ (-\infty)\f$) and round up 
(or toward \f$ (+\infty)\f$). When not specified explicitly, the 
operations use the default rounding mode, which is in practice a 
variable local to each execution thread. The default rounding mode 
can be set to any of the four rounding modes (initially, it is set 
to nearest). To specify rounding modes for operations, the type 
used is `std::float_round_style`. 

This type is `ImplicitInteroperable` with `Gmpz`, <TT>long</TT>, 
<TT>unsigned long</TT>, <TT>int</TT>, <TT>double</TT> and <TT>long double</TT>. 

\cgalModels `FieldWithKthRoot` 
\cgalModels `RealEmbeddable` 

\cgalHeading{Comparisons}

Comparison operators `==`, `!=`, `>`, `<`, 
`>=` and `<=` are also overloaded. A `Gmpfr` can 
be compared with other `Gmpfr`, as well as with a `Gmpz`, 
`long`, `unsigned long`, `int`, `double` or 
`long double`. It is worth noting that the numbers are never 
converted nor rounded before comparison. In the case where one of 
the compared numbers is `NaN`, the `erange` flag is set. 

\cgalHeading{Implementation}

Since the <span class="textsc">Mpfr</span> library can be compiled to be thread-safe, this interface 
is designed to keep the thread-safety. 

`Gmpfr`s are reference counted. This behavior may be changed, by 
setting the flag <TT>CGAL_GMPFR_NO_REFCOUNT</TT>. A non-reference-counted 
class is slightly more efficient in case the implementation does not need 
to copy numbers (this is not usually the case). Nevertheless, setting this 
flag may be useful for debugging purposes. 

\sa `RealEmbeddable` 
\sa `FieldWithKthRoot` 

*/
class Gmpfr {
public:

/// \name Types 
/// @{

/*!
Type representing the precision (number of bits 
used to represent the mantissa) of a number. 
*/ 
typedef unspecified_type Precision_type; 

/// @} 

/// \name Creation 
/// Note that all constructors can be called with two optional
/// parameters. One can specify as second parameter the rounding mode
/// desired for the conversion from the source number and as a third
/// parameter the precision with which this `Gmpfr` will be
/// created. If only one optional parameter is specified, it can be
/// either the rounding mode or the precision. If no optional
/// parameters are specified, the precision of the created object is
/// chosen in such a way that the conversion is exact (i.e., no
/// rounding is performed). 
///
/// These optional parameters, along with
/// other functions which will be explained below, allow users to
/// control the rounding and precision. For example, being <TT>z</TT>
/// a `Gmpz`, <TT>Gmpfr g(z,53,std::round_toward_neg_infinity)</TT>
/// will construct a `Gmpfr` <TT>g</TT> having as value the biggest
/// 53-bit floating-point number that is equal or smaller than to
/// <TT>z</TT>.
/// @{

/*!
Creates an uninitialized `Gmpfr` `f`. 
*/ 
Gmpfr(); 

/*!
Copy constructor. The copied object inherits the precision of 
`n`, and thus it is not rounded. 
*/ 
Gmpfr(const Gmpfr& n); 

/*!
Creates a `Gmpfr`, initialized with the value of `si`. 
*/ 
Gmpfr(long si); 

/*!
Creates a `Gmpfr`, initialized with the value of `ui`. 
*/ 
Gmpfr(unsigned long ui); 

/*!
Creates a `Gmpfr`, initialized with the value of `i`. 
*/ 
Gmpfr(int i); 

/*!
Creates a `Gmpfr`, initialized with the value of `d`. 
*/ 
Gmpfr(double d); 

/*!
Creates a `Gmpfr`, initialized with the value of `ld`. 
*/ 
Gmpfr(long double ld); 

/*!
Creates a `Gmpfr`, initialized with the value of `z`. 
*/ 
Gmpfr(const Gmpz &z); 

/*!
Creates a `Gmpfr`, initialized with the value of `zf`. 
*/ 
Gmpfr(const Gmpzf &zf); 

/*!
Creates a `Gmpfr`, initialized with the value of 
`ie.first`\f$ \times2^{\mathrm{ie.second}} \f$ . 
*/ 
Gmpfr(const std::pair<Gmpz,long> &ie); 

/*!
This returns the current precision used in `Gmpfr` 
creation by default. 
*/ 
static Precision_type get_default_precision(); 

/*!
This function sets the default <span class="textsc">Mpfr</span> precision to p, and returns 
the old one. 
*/ 
static Precision_type set_default_precision(Precision_type p); 

/// @} 

/*! \name Controlling the Precision 

Each Gmpfr object has a precision associated to it. The precision is
the amount of bits needed to represent the mantissa. <span
class="textsc">Mpfr</span> has a default precision value, which can be
controlled by static functions of the Gmpfr class (in practice, this
default value is a variable local to each execution thread). There are
also functions to get and set the precision of each Gmpfr object.
*/

/// @{

/*!
Returns the precision of `f`. 
*/ 
Precision_type get_precision()const; 

/*!
Returns the value of `f`, rounded with precision `p` 
in the direction `r`. 
*/ 
Gmpfr round(Precision_type p, std::float_round_style r)const; 

/*!
This function returns the current rounding mode used by <span class="textsc">Mpfr</span>. 
*/ 
static std::float_round_style get_default_rndmode(); 

/*!
This function sets the <span class="textsc">Mpfr</span> rounding mode to `r` and returns 
the old one. 
*/ 
static std::float_round_style set_default_rndmode(std::float_round_style r); 

  /// @}

/*! \name Flags

\sc{Mpfr} provides some flags to know whether
performed operations were exact or not, or they incurred in overflow
or underflow, if the exponent is out of range, or the result was `NaN`
(not-a-number). One can clear the flags before a set of operations and
inspect them afterward, in order to see if something unexpected
happened during the operations. The static functions used to handle
flags are:
*/

/// @{
/*!
Clears all the flags set by <span class="textsc">Mpfr</span>(they are not cleared 
automatically). 
*/ 
static void clear_flags(); 

/*!
Shows whether an operation incurred in underflow. 
*/ 
static bool underflow_flag(); 

/*!
Shows whether an operation incurred in overflow. 
*/ 
static bool overflow_flag(); 

/*!
Shows whether the result of an operation was `NaN`. 
*/ 
static bool nan_flag(); 

/*!
Shows whether an operation was inexact. 
*/ 
static bool inex_flag(); 

/*!
Returns `true` iff a range error occurred. Such an exception 
occurs when some function which does not return a `Gmpfr` 
has an invalid result. For example, this flag will be set if 
one of the operands of a comparison is `NaN`. 
*/ 
static bool erange_flag(); 

/// @}

/*! \name Arithmetic Operations

 Arithmetic operators `+` , `-` , and `/` are overloaded, but special
 care must be taken when applying them. The precision of an operation
 between two `Gmpfr`s is defined as the maximum of the operands
 precision and the default precision. The second operand of the former
 operations can be a `Gmpfr`, `int`, `long`, `unsigned`, `unsigned
 long`, or `Gmpz`. The precision of an operation between a `Gmpfr` and
 a number of another type is defined as the maximum between the
 number's precision and the default precision. To specify the rounding
 mode and/or the precision to perform an operation, this class
 provides the four static functions `add`, `sub`, `mul` and
 `div`. Only one of them is shown here, since their interfaces are
 similar: When the precision is not specified in this family of
 functions, it is defined as in the overloaded operators. When the
 rounding mode is not specified, the default is used.
*/

/// @{

/*!
*/ 
static Gmpfr add (const Gmpfr &a,const Gmpfr &b); 

/*!
*/ 
static Gmpfr add (const Gmpfr &a,const Gmpfr &b, std::float_round_style r); 

/*!
*/ 
static Gmpfr add (const Gmpfr &a,const Gmpfr &b,Precision_type p); 

/*!
*/ 
static Gmpfr add (const Gmpfr &a,const Gmpfr &b,Precision_type p, std::float_round_style r); 

/*!
Returns the absolute value of `f`, rounded with precision 
`p` in the direction `r`. If `p` is not specified, 
the precision used is the maximum between `f`'s precision 
and the default. 
*/ 
Gmpfr abs(Precision_type p, std::float_round_style r=get_default_rndmode())const; 

/*!
Returns the square root of `f`, rounded with precision 
`p` in the direction `r`. If `p` is not specified, 
the precision used is the maximum between `f`'s precision 
and the default. 
*/ 
Gmpfr sqrt(Precision_type p, std::float_round_style r=get_default_rndmode())const; 

/*!
Returns the k-th root of `f`, rounded with precision 
`p` in the direction `r`. If `p` is not specified, 
the precision used is the maximum between `f`'s precision 
and the default. 
*/ 
Gmpfr kthroot(int k, Precision_type p, std::float_round_style r=get_default_rndmode())const; 

/*!
Returns the square of `f`, rounded with precision `p` in 
the direction `r`. If `p` is not specified, the precision 
used is the maximum between `f`'s precision and the default. 
*/ 
Gmpfr square(Precision_type p, std::float_round_style r=get_default_rndmode())const; 

/*!
Returns a double precision approximation of `f` using the 
rounding mode `r`. 
*/ 
double to_double(std::float_round_style r=get_default_rndmode()); 

/*!
Returns an interval of doubles which contains `f`. If a 
rounded endpoint does not fit in a double, the double is set to plus 
or minus infinity and the `overflow` or `underflow` flag. 
*/ 
std::pair<double,double> to_interval(); 

/*!
Returns the pair \f$ (d,e) \f$ such that \f$ 0.5 \le|d| < 1 \f$ and 
\f$ d \times2^e \f$ equals `f` rounded to double precision, 
using the rounding mode `r`. If `f` is `NaN` or 
infinity, then the corresponding double is returned, leaving 
the exponent undefined and setting the appropriate error flag. 
*/ 
std::pair<double,long> to_double_exp 
(std::float_round_style r=get_default_rndmode()); 

/*!
Returns \f$ ((m,M),e) \f$ such that \f$ m \times2^e \le f \le M \times2^e \f$. If `f` is `NaN` or infinity, then 
the corresponding doubles are returned, leaving the exponent 
undefined and setting the appropriate error flag. 
*/ 
std::pair<std::pair<double,double>,long> to_interval_exp(); 

/*!
Returns a pair of integers \f$ (m,e) \f$, such that 
\f$ f = m \times2^e \f$. Note that the returned value of \f$ m\f$ 
is not necessarily the smallest possible value of \f$ m\f$ (that is, 
it might be that \f$ 2|m\f$). 
*/ 
std::pair<Gmpz,long> to_integer_exp(); 

/// @} 

/// \name Query Functions 
/// @{

/*!
Returns the sign of `f`. 
*/ 
Sign sign(); 

/*!
Returns `true` iff `f` is zero. 
*/ 
bool is_zero(); 

/*!
Returns `true` iff `f` is one. 
*/ 
bool is_one(); 

/*!
Returns `true` iff `f` is NaN (not-a-number). 
*/ 
bool is_nan(); 

/*!
Returns `true` iff `f` is plus or minus infinity. 
*/ 
bool is_inf(); 

/*!
Returns `true` iff `f` is a valid number. 
*/ 
bool is_number(); 

/*!
Returns `true` iff `f` is the square of a number 
representable by an object of this type. 
*/ 
bool is_square(); 

/*!
Returns `true` iff `f` is the square of a number 
representable by an object of this type, computing and storing it 
in `y`. 
*/ 
bool is_square(const Gmpfr &y); 


/// @}

}; /* end Gmpfr */

/*!
Reads a floating-point number from `in`. The number 
\f$ M \times2^E\f$ must be in the form \f$ MeE\f$, where the mantissa 
\f$ M\f$ and the exponent \f$ E\f$ are integers in base 10. 
\relates Gmpfr 
*/ 
std::istream& operator>>(std::istream& in, Gmpfr& f); 

/*!
If the ostream `out` is in pretty-print mode, writes a decimal 
approximation of `f` to `out`. Otherwise, writes `f` to 
`out` in the form \f$ MeE\f$, where \f$ M\f$ is its mantissa and 
\f$ E\f$ is its exponent, both in base 10. 
\relates Gmpfr 
*/ 
std::ostream& operator<<(std::ostream& out, const Gmpfr& f); 

} /* end namespace CGAL */
