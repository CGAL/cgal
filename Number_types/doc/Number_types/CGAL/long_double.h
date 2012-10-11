
/*!
\ingroup nt_builtin

The fundamental type `long double` is a model of the concepts `RealEmbeddable` 
and `FieldWithSqrt`. Due to rounding errors and overflow `long double` is 
considered as not exact. 

\models ::FieldWithSqrt 
\models ::RealEmbeddable 

*/

class long_double {
}; /* end long double */

/*! 
Determines whether the argument represents a value in \f$ \R\f$. 
\relates long_double 
*/ 
bool CGAL::is_finite(long double x); 

