namespace CGAL {

/*!
\ingroup nt_rrational

computes the rational number that equals `d`.

The function `to_rational` computes the rational number representing a 
given double precision floating point number. 

\sa `CGAL::simplest_rational_in_interval<Rational>(double d1, double d2)`
*/
template <typename Rational>
Rational to_rational(double d);

} /* namespace CGAL */

