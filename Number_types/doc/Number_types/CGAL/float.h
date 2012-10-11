

/*!
\ingroup nt_builtin

The fundamental type `float` is a model of the concepts `RealEmbeddable` and 
`FieldWithSqrt`. Due to rounding errors and overflow `float` is 
considered as not exact. 

\models ::FieldWithSqrt 
\models ::RealEmbeddable 

*/
class float {

}; /* end float */

namespace CGAL {
/*! 
Determines whether the argument represents a value in \f$ \R\f$. 
\relates float 
*/ 
bool CGAL::is_finite(float x); 

} /* end namespace CGAL */
