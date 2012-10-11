

/*!
\ingroup nt_builtin

The fundamental type `double` is a model of the concepts `RealEmbeddable`
and `Field`. Due to rounding errors and overflow `double` is 
considered as not exact. 

\models ::FieldWithSqrt 
\models ::RealEmbeddable 

*/
class double {
}; /* end double */

namespace CGAL {
/*! 
Determines whether the argument represents a value in \f$ \R\f$. 
\relates double 
*/ 
bool CGAL::is_finite(double x); 

} // end namespace

