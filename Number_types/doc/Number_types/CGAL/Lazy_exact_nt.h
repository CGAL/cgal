
namespace CGAL {

/*!
\ingroup nt_cgal

An object of the class `Lazy_exact_nt<NT>` is able to represent any
real embeddable number which `NT` is able to represent.
The idea is that `Lazy_exact_nt<NT>` works exactly like `NT`, except
that it is expected to be faster because it tries to only compute an
approximation of the value, and only refers to `NT` when needed.
The goal is to speed up exact computations done by any exact but slow
number type `NT`.

The function `to_double()` can be used to get a double approximation
of the represented number. Note that two subsequent calls to this
function on the same number of type `Lazy_exact_nt<NT>` might not return
the same value as the exact representation might have been computed
between the two calls, thus refining the double approximation. If you
want to avoid this behavior, you need to first call `exact()`
(losing the benefit of the laziness if done systematically).

\tparam NT must be a model of concept `RealEmbeddable`, and at
least model of concept `IntegralDomainWithoutDivision`.

Note that some filtering mechanism is available at the predicate level
using `Filtered_predicate` and `Filtered_kernel`.

\cgalModels{IntegralDomainWithoutDivision same as `NT`,RealEmbeddable,
             Fraction, if `NT` is a `Fraction`}

\cgalHeading{Example}

\code
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Quotient.h>

typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > NT;
typedef CGAL::Cartesian<NT> K;
\endcode

*/
template< typename NT >
class Lazy_exact_nt {
public:

/// \name Creation
/// @{

/*!
introduces an uninitialized variable `m`.
*/
Lazy_exact_nt();

/*!
introduces the value \a x, of any built-in arithmetic type (`int`, `double`, etc) (works only if `NT` has a constructor from this type too).
*/
Lazy_exact_nt(BuiltIn i);

/*!
introduces the value `n`.
*/
Lazy_exact_nt(NT n);

/*!
introduces the value `n`. `NT1` needs to be convertible to `NT`
(and this conversion will only be done if necessary).
*/
template <class NT1> Lazy_exact_nt(Lazy_exact_nt<NT1> n);

/// @}

/// \name Operations
/// @{

/*!
returns the corresponding NT value.
*/
NT exact();

/*!
returns an interval containing the
exact value.
*/
Interval_nt<false> approx();

/*!
returns an interval containing the
exact value.
*/
Interval_nt<true> interval();

/*!
specifies the relative precision that `to_double()` has to fulfill.
The relative precision is thread local, and the default value is \f$ 10^{-5}\f$.

\pre `d>0` and `d<1`.

*/
static void set_relative_precision_of_to_double(double d);

/*!
returns the relative precision that `to_double()` currently fulfills.
*/
static double get_relative_precision_of_to_double();

/// @}

}; /* end Lazy_exact_nt */

/*!
writes `m` to ostream `out` in an interval format.
\relates Lazy_exact_nt
*/
std::ostream& operator<<(std::ostream& out, const Lazy_exact_nt<NT>& m);

/*!
reads a `NT` from `in`, then converts it to a `Lazy_exact_nt<NT>`.
\relates Lazy_exact_nt
*/
std::istream& operator>>(std::istream& in, Lazy_exact_nt<NT>& m);


} /* end namespace CGAL */
