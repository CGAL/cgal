
namespace CGAL {

/*!
\ingroup nt_cgal

The class `Interval_nt` provides an interval arithmetic number type.

This section describes briefly what interval arithmetic is, its implementation
in \cgal, and its possible use by geometric programs.
The main reason for having interval arithmetic in \cgal is its integration
into the filtered robust and fast predicates scheme, but we also provide a
number type so that you can use it separately if you find any use for it,
such as interval analysis, or to represent data with tolerance...

The purpose of interval arithmetic is to provide an efficient way to bound
the roundoff errors made by floating point computations.
You can choose the behavior of your program depending on these errors.
You can find more theoretical information on this topic in
\cgalCite{cgal:bbp-iayed-01}.

Interval arithmetic is a large concept and we will only consider here a
simple arithmetic based on intervals whose bounds are <I>double</I>s.
So each variable is an interval representing any value inside the interval.
All arithmetic operations (+, -, \f$ *\f$, \f$ /\f$, \f$ \sqrt{}\f$, `square()`,
`min()`, `max()` and `abs()`) on intervals preserve the inclusion.
This property can be expressed by the following formula (\f$ x\f$ and \f$ y\f$ are
real, \f$ X\f$ and \f$ Y\f$ are intervals, \f$ \mathcal{OP}\f$ is an arithmetic operation):

\f[
\forall\ x \in X, \forall\ y \in Y, (x\ \mathcal{OP}\ y)
\in (X\ \mathcal{OP}\ Y)
\f]

For example, if the final result of a sequence of arithmetic operations is
an interval that does not contain zero, then you can safely determine its sign.

\cgalHeading{Parameters}

The template parameter `Protected` is a Boolean parameter, which defaults
to `true`. It provides a way to select faster computations by avoiding
rounding mode switches, at the expense of more care to be taken by the user
(see below). The default value, `true`, is the safe way, and takes care of
proper rounding mode changes. When specifying `false`, the user has to
take care about setting the rounding mode towards plus infinity before
doing any computations with the interval class. He can do so using the
`Protect_FPU_rounding` class for example.

\cgalModels `FieldWithSqrt`
\cgalModels `RealEmbeddable`

\cgalHeading{Example}

Protecting an area of code that uses operations on the class
`Interval_nt_advanced` can be done in the following way:

\code
{
Interval_nt_advanced::Protector P;
... // The code to be protected.
}
\endcode

The basic idea is to use the directed rounding modes specified by the
<I>IEEE 754</I> standard, which are implemented by almost all processors
nowadays.
It states that you have the possibility, concerning the basic floating point
operations (\f$ +,-,*,/,\sqrt{}\f$) to specify the rounding mode of each operation
instead of using the default, which is set to 'round to the nearest'.
This feature allows us to compute easily on intervals. For example, to
add the two intervals [a.i;a.s] and [b.i;b.s], compute \f$ c.i=a.i+b.i\f$ rounded
towards minus infinity, and \f$ c.s=a.s+b.s\f$ rounded towards plus infinity, and
the result is the interval [c.i;c.s]. This method can be extended easily to
the other operations.

The problem is that we have to change the rounding mode very often, and the
functions of the C library doing this operation are slow and not portable.
That's why assembly versions are used as often as possible.
Another trick is to store the opposite of the lower bound, instead of the
lower bound itself, which allows us to never change the rounding mode inside
simple operations. Therefore, all basic operations, which are in the class
`Interval_nt_advanced` assume that the rounding mode is set to
'round to infinity', and everything works with this correctly set.

So, if the user needs the speed of `Interval_nt_advanced`, he must
take care of setting the rounding mode to 'round to infinity' before each
block of operations on this number type. And if other operations might be
affected by this, he must take care to reset it to 'round to the nearest'
before they are executed.

Notes:

<UL>
<LI>On Intel platforms (with any operating system and compiler), due to a
misfeature of the floating point unit, which does not handle exactly
IEEE compliant operations on doubles, we are forced to use a workaround
which slows down the code, but is only useful when the intervals can
overflow or underflow. If you know that the intervals will never
overflow nor underflow for your code, then you can disable this
workaround with the flag `CGAL_IA_NO_X86_OVER_UNDER_FLOW_PROTECT`.
Other platforms are not affected by this flag.
<LI>When optimizing, compilers usually propagate the value of variables when
they know it's a constant. This can break the interval routines because
the compiler then does some floating point operations on these constants
with the default rounding mode, which is wrong. This kind of problem
is avoided by stopping constant propagation in the interval routines.
However, this solution slows down the code and is rarely useful, so you
can disable it by setting the flag
`CGAL_IA_DONT_STOP_CONSTANT_PROPAGATION`.
</UL>

*/
template< typename Protected >
class Interval_nt {
public:

/// \name Types
/// The class `Interval_nt` defines the following types:
/// @{

/*!
The type of the bounds of the interval.
*/
typedef double value_type;

/*!
The type of the
exceptions raised when uncertain comparisons are performed.
*/
typedef Uncertain_conversion_exception unsafe_comparison;

/*!
A type whose default constructor and destructor allow
to protect a block of code from FPU rounding modes necessary for the
computations with `Interval_nt<false>`. It does nothing for
`Interval_nt<true>`. It is implemented as `Protect_FPU_rounding<!Protected>`.
*/
typedef unspecified_type Protector;

/// @}

/// \name Creation
/// @{

/*!
introduces a small interval containing \a i (possibly a point).
*/
Interval_nt(long long i);

/*!
introduces the interval [`d`;`d`].
*/
Interval_nt(double d);

/*!
introduces the interval [`i`;`s`].
*/
Interval_nt(double i, double s);

/*!
introduces the interval [`p.first`;`p.second`].
*/
Interval_nt(std::pair<double, double> p);

/// @}

/// \name Operations
/// All functions required by a class to be considered as a \cgal number type (see \ref Numbertype) are present, as well as the utility functions, sometimes with a particular semantic which is described below. There are also a few additional functions.
/// @{

/*!
returns [\f$ -\infty\f$;\f$ +\infty\f$] when the denominator contains 0.
*/
Interval_nt operator/(Interval_nt J);


/*!
returns the lower bound of the interval.
*/
double inf();

/*!
returns the upper bound of the interval.
*/
double sup();

/*!
returns whether both bounds are equal.
*/
bool is_point();

/*!
returns whether both intervals have
the same bounds.
*/
bool is_same(Interval_nt J);

/*!
returns whether both intervals
have a non empty intersection.
*/
bool do_overlap(Interval_nt J);

/// @}

/*! \name Implementation

The operations on `Interval_nt` with the default parameter `true`, are
automatically protected against rounding modes, and are thus slower
than those on `Interval_nt_advanced`, but easier to use. Users that
need performance are encouraged to use `Interval_nt_advanced`
instead. Changing the rounding mode affects all floating point
computations, and might cause problems with parts of your code, or
external libraries (even \cgal), that expect the rounding mode to be
the default (round to the nearest).

We provide two interfaces to
change the rounding mode. The first one is to use a protector object
whose default constructor and destructor will take care of changing
the rounding mode. The protector is implemented using
`Protect_FPU_rounding`. The second one is the following set of
functions. The macros `CGAL_FE_TONEAREST`, `CGAL_FE_TOWARDZERO`,
`CGAL_FE_UPWARD` and `CGAL_FE_DOWNWARD` are the values corresponding
to the rounding modes.

*/

/// @{

/*!
The type used by the following functions to deal with rounding modes.
This is usually an `int`.
*/
typedef int FPU_CW_t;

/// @}

}; /* end Interval_nt */

/*!
returns [0;\f$ \sqrt{upper\_bound(I)}\f$] when only the lower bound is negative (expectable
case with roundoff errors), and is unspecified when the upper bound also is
negative (unexpected case).
\relates Interval_nt
*/
Interval_nt sqrt(Interval_nt I);

/*!
returns the
middle of the interval, as a double approximation of the interval.
\relates Interval_nt
*/
double to_double(Interval_nt I);

  /*! \name Comparisons
The comparison operators (\f$ <\f$, \f$ >\f$, \f$ <=\f$, \f$ >=\f$, \f$ ==\f$, \f$ !=\f$, `sign()` and `compare()`) have the following semantic: it is the intuitive one when for all couples of values in both intervals, the comparison is identical (case of non-overlapping intervals). This can be expressed by the following formula (\f$ x\f$ and \f$ y\f$ are real, \f$ X\f$ and \f$ Y\f$ are intervals, \f$ \mathcal{OP}\f$ is a comparison operator): \f[ \left(\forall x \in X, \forall y \in Y, (x\ \mathcal{OP}\ y) = true\right) \Rightarrow (X\ \mathcal{OP}\ Y) = true \f] and \f[ \left(\forall x \in X, \forall y \in Y, (x\ \mathcal{OP}\ y) = false\right) \Rightarrow (X\ \mathcal{OP}\ Y) =false \f] Otherwise, the comparison is not safe, and we specify this by returning a type encoding this uncertainty, namely using `Uncertain<bool>` or `Uncertain<Sign>`, which can be probed for uncertainty explicitly, and which has a conversion to the normal type (e.g. `bool`) which throws an exception when the conversion is not certain. Note that each failed conversion increments a profiling counter (see `CGAL_PROFILE`), and then throws the exception of type `unsafe_comparison`.
   */

 /// @{

/*!

\relates Interval_nt
*/
Uncertain<bool> operator<(Interval_nt i, Interval_nt j);

/*!

\relates Interval_nt
*/
Uncertain<bool> operator>(Interval_nt i, Interval_nt j);

/*!

\relates Interval_nt
*/
Uncertain<bool> operator<=(Interval_nt i, Interval_nt j);

/*!

\relates Interval_nt
*/
Uncertain<bool> operator>=(Interval_nt i, Interval_nt j);

/*!
\relates Interval_nt
*/
Uncertain<bool> operator==(Interval_nt i, Interval_nt j);


/*!
\relates Interval_nt
*/
Uncertain<bool> operator!=(Interval_nt i, Interval_nt j);


/*!
\relates Interval_nt
*/
Uncertain<Comparison_result>
compare(Interval_nt i, Interval_nt j);

 /// @}

/*!
\relates Interval_nt
*/
Uncertain<Sign> sign(Interval_nt i);

/*!
sets the rounding mode to `R`.
\relates Interval_nt
*/
void FPU_set_cw (FPU_CW_t R);

/*!
returns the current rounding mode.
\relates Interval_nt
*/
FPU_CW_t FPU_get_cw (void);

/*!
sets the rounding mode to `R` and returns the old one.
\relates Interval_nt
*/
FPU_CW_t FPU_get_and_set_cw (FPU_CW_t R);


/*!
\ingroup nt_cgal
This typedef (at namespace \cgal scope) exists for backward compatibility,
as well as removing the need to remember the Boolean value for the template
parameter.
*/
typedef Interval_nt<false> Interval_nt_advanced;

} /* end namespace CGAL */
