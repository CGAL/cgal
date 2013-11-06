
namespace CGAL {

/*!
\ingroup nt_gmp

An object of the class `Gmpfi` is a closed interval, with endpoints 
represented as `Gmpfr` floating-point numbers. An interval can have 
finite or infinite endpoints and its meaning is straightforward. It can 
also have one (or both) `NaN` endpoint(s): this indicates that an 
invalid operation has been performed and that the resulting interval has no 
mathematical meaning. 

All the operations of `Gmpfi` were designed in such a way that the 
mathematical correct result is always contained in the resulting interval. 

This type is `ImplicitInteroperable` with `Gmpfr`, `Gmpz`, 
`Gmpq`, <TT>long</TT>, <TT>unsigned long</TT>, <TT>int</TT>, <TT>double</TT> 
and <TT>long double</TT>. 

\cgalModels `FieldWithKthRoot` 
\cgalModels `RealEmbeddable` 

\cgalHeading{Implementation}

All interval operations are performed by the <span
class="textsc">Mpfi</span> library. The class `Gmpfi` is not reference
counted, but its members are.

The default precision of `Gmpfi` is local to each thread and independent of
the default precision of `Gmpfr` (in contrast to the behaviour of the <span
class="textsc">Mpfi</span> and <span class="textsc">Mpfr</span> libraries,
which share a default precision).

\sa `CGAL::Gmpfr` 
\sa `CGAL::Interval_nt` 
\sa `CGAL::Uncertain` 
\sa `RealEmbeddable` 
\sa `FieldWithKthRoot` 

*/
class Gmpfi {
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
/// All constructors accept as an optional last argument a precision,
/// which can be used to specify the precision of
/// the `Gmpfr` endpoints. If none is specified, the default precision
/// will be used. As the endpoints are represented with a fixed number
/// of bits, they may need to be rounded. In this case, the number
/// from which the `Gmpfi` was constructed is guaranteed to be
/// included in the constructed interval.
/// @{

/*!
creates an uninitialized `Gmpfi` interval `i`. 
*/ 
Gmpfi(Precision_type p=get_default_precision()); 

/*!
creates a `Gmpfi` initialized with the value of `t`. 
`T` is `Gmpfr`, `Gmpq`, or any type from which 
`Gmpfr` can be constructed from. The rounding of the 
endpoints will guarantee that `t` is included in `i`. 
*/ 
template <class T> 
Gmpfi(const T& t,Precision_type p=get_default_precision()); 

/*!
creates a `Gmpfi` initialized with endpoints `left` 
and `right`. The rounding of the endpoints will guarantee 
that ([`left`,`right`]) is included in `i`. 
*/ 
Gmpfi(const Gmpfr &left, 
const Gmpfr &right, 
Precision_type p=get_default_precision()); 

/*!
creates a `Gmpfi` initialized with endpoints 
`endpoints.first` and `endpoints.second`. `L` and 
`R` are types from which `Gmpfr` can be constructed 
from. The rounding of the endpoints will guarantee that 
([`endpoints.first`,`endpoints.second`]) is included in 
`i`. 
*/ 
template<class L, class R> 
Gmpfi(const std::pair<L,R> &endpoints, 
Precision_type p=get_default_precision()); 

/// @} 

/// \name Operations 

/// @{

/*!
Returns the smallest (or <I>inferior</I>) `Gmpfr` endpoint of 
`i`. 
*/ 
Gmpfr inf() const; 

/*!
Returns the largest (or <I>superior</I>) `Gmpfr` endpoint of 
`i`. 
*/ 
Gmpfr sup()const; 

/*!
Returns the default precision. 
*/ 
static Precision_type get_default_precision(); 

/*!
Sets the default precision to `prec` and returns the 
old value. 
*/ 
static Precision_type set_default_precision(Precision_type prec); 

/*!
Returns the precision of `i`. 
*/ 
Precision_type get_precision()const; 

/*!
Returns the value of the number, rounded with precision `p`. 
*/ 
Gmpfi round(Precision_type p)const; 

/// @} 

/*! \name Arithmetic Operations 

 
Arithmetic operators `+`, `-`, and `/` are overloaded, but special
care must be taken when applying them.  The precision of an operation
between two `Gmpfi`s is defined as the maximum of the operands
precision and the default precision.

The second operand of the former operations can be a `Gmpfi`, `Gmpfr`,
`int`, `long`, `unsigned`, `unsigned long`, `Gmpz` or `Gmpq`. The
precision of an operation between a `Gmpfi` and a number of another
type is defined as the `Gmpfi`'s precision (even when operating with a
`Gmpfr`).

To specify the rounding mode and/or the precision to perform an
operation, this class provides the four static functions `add`, `sub`,
`mul` and `div`. Only one of them is shown here, since their
interfaces are similar: When the precision is not specified in this
family of functions, it is defined as in the overloaded operators.
*/
 
/// @{
static Gmpfi add (const Gmpfi &a,const Gmpfi &b,Precision_type p=0); 

/*!
Returns the absolute value of `i`, with precision `p`. 
If `p` is not specified, the precision used is the maximum 
between `i`'s precision and the default. 
*/ 
Gmpfi abs(Precision_type p)const; 

/*!
Returns the square root of `i`, with precision `p`. 
If `p` is not specified, the precision used is the maximum 
between `i`'s precision and the default. 
*/ 
Gmpfi sqrt(Precision_type p)const; 

/*!
Returns the k-th root of `i`, with precision `p`. 
If `p` is not specified, the precision used is the maximum 
between `i`'s precision and the default. 
*/ 
Gmpfi kthroot(int k,Precision_type p)const; 

/*!
Returns the square of `i`, with precision `p`. If 
`p` is not specified, the precision used is the maximum 
between `i`'s precision and the default. 
*/ 
Gmpfi square(Precision_type p)const; 

/*!
Returns an interval of doubles which contains `i`. If a 
rounded endpoint does not fit in a double, sets its value to plus 
or minus infinity and the `overflow` or `underflow` flag. 
*/ 
std::pair<double,double> to_interval()const; 

/*!
Returns \f$ (m,e)\f$ such that \f$ m \times2^e\f$ is the center of 
`i`, rounded to nearest. If one of the endpoints of `i` is 
`NaN` or infinity, then the corresponding double is returned, 
leaving the exponent undefined and setting the appropriate 
error flag. 
*/ 
std::pair<double,long> to_double_exp()const; 

/*!
Returns \f$ ((m_1,m_2),e)\f$, such that \f$ [m_1 \times2^e,m_2 
\times2^e]\f$ contains `i`. If one of the endpoints of 
`i` is `NaN` or infinity, then the corresponding doubles 
are returned, leaving the exponent undefined and setting the 
appropriate error flag. 
*/ 
std::pair<std::pair<double,double>,long> to_interval_exp()const; 

/// @} 
/*! \name Comparisons 

The semantics of the comparison operators is the same than on `Interval_nt<Protected>`. The result of the comparison is always an `Uncertain<bool>` (this type is convertible to `bool`, but may throw an exception). If compared intervals have no common points, the result is `true` or `false`; otherwise, `Uncertain::indeterminate()` will be returned. 

In the same way, we can explain the semantics of `Uncertain<Comparison_result>` and `Uncertain<Sign>`. With the semantics described above, this class provides comparisons between `Gmpfi` and `Gmpfi`, `Gmpfr`, `long`, `unsigned long`, `int`, `double`, `Gmpz` and `Gmpq`. Comparison operators `==`, `!=`, `>`, `<`, `>=` and `<=` are overloaded. 

The class provides also functions to test efficiently some special kinds of comparisons:
*/

/// @{

/*!
Returns `true` iff left endpoints of `i` and 
`j` are equal and right endpoints of them are also equal. Note 
that this does not mean equality between `i` and `j`. 
*/ 
bool is_same(const Gmpfi &j)const; 

/*!
Returns `true` iff `i` and `j` overlap, 
i.e., iff they have points in common. 
*/ 
bool do_overlap(const Gmpfi &j)const; 

/*!
If `i` and `j` do not overlap, this function returns 
the result of the comparison. Otherwise, it returns 
`indeterminate`. 
*/ 
Uncertain<Comparison_result> compare(const Gmpfi &j)const; 

/// @} 

/// \name Query Functions 
/// @{

/*!
Returns `true` iff both endpoints are equal. 
*/ 
bool is_point()const; 

/*!
Returns `true` iff at least one of the endpoints is 
`NaN`. 
*/ 
bool is_nan()const; 

/*!
Returns `true` iff at least one of the endpoints is plus or 
minus infinity. 
*/ 
bool is_inf()const; 

/*!
Returns `true` iff `i` is a bounded interval, i.e., its 
endpoints are neither invalid nor infinite. 
*/ 
bool is_number()const; 

/*!
Returns `true` if both endpoints are zero, `false` if 
the interval does not contain zero and `indeterminate` 
otherwise. 
*/ 
Uncertain<bool> is_zero()const; 

/*!
Returns `true` if both endpoints are one, `false` 
if the interval does not contain one and `indeterminate` 
otherwise. 
*/ 
Uncertain<bool> is_one()const; 

/*!
If all numbers contained in the interval have the same sign, this 
function returns it. Otherwise it returns `indeterminate`. 
*/ 
Uncertain<Sign> sign()const; 

/*!
Returns `true` if all numbers contained in the interval 
are positive, false if all of them are negative or zero and 
`indeterminate` otherwise. 
*/ 
Uncertain<bool> is_positive()const; 

/*!
Returns `true` if all numbers contained in the interval 
are negative, false if all of them are positive or zero and 
`indeterminate` otherwise. 
*/ 
Uncertain<bool> is_negative()const; 


/// @}

}; /* end Gmpfi */

/*!
Reads `i` from `is`. `is` must have the form 
`[inf,sup]`, where `inf` and `sup` have valid 
`Gmpfr` input formats. 
\relates Gmpfi 
*/ 
std::istream& operator>>(std::istream &is,Gmpfi i); 

/*!
Writes `i` to `os`, in the form `[i.inf(),i.sup()]`. 
The endpoints are written according to the `Gmpfr` formatting. 
\relates Gmpfi 
*/ 
std::ostream& operator<<(std::ostream &os,const Gmpfi &i); 
} /* end namespace CGAL */
