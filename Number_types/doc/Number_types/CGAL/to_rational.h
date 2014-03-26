namespace CGAL {

/*!
\ingroup nt_rrational

computes the rational number that equals `d`.

Computes the rational number representing a 
given double precision floating point number. 

\sa `CGAL::simplest_rational_in_interval()`
*/
template <typename Rational>
Rational to_rational(double d);

} /* namespace CGAL */

