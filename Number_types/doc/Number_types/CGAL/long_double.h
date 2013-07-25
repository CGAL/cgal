
/*!
\file long_double.h
\ingroup nt_builtin

This header provides all necessary functions so the fundamental type
`long double` is a model of the concepts `RealEmbeddable` and
`FieldWithSqrt`. Due to rounding errors and overflow `long double` is
considered as not exact.
*/

namespace CGAL {

/*!
Determines whether the argument represents a value in \f$ \mathbb{R}\f$. 
*/ 
bool is_finite(long double x); 

} // CGAL


