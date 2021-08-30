namespace CGAL {

/*!
\ingroup nt_rrational

computes the rational number with the smallest denominator in the
interval `[d1,d2]`.

Computes the simplest rational number in an
interval of two `double` values.

\cgalHeading{Implementation}

See Knuth, "Seminumerical algorithms", page 654, answer to exercise
4.53-39.

\sa `CGAL::to_rational()`

*/
template <typename Rational>
Rational simplest_rational_in_interval(double d1, double d2);

} /* namespace CGAL */

