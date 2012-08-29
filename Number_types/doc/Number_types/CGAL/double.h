
namespace CGAL {

/*!
\ingroup nt_builtin

The fundamental type `double` is an `RealEmbeddable` 
`Field`. Due to rounding errors and overflow `double` is 
considered as not exact. 

\models ::FieldWithSqrt 
\models ::RealEmbeddable 

*/
class double {
}; /* end double */
/*! 
Determines whether the argument represents a value in \f$ \R\f$. 
\relates double 
*/ 
bool CGAL::is_finite(double x); 


} /* end namespace CGAL */
