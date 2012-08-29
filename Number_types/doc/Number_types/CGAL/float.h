namespace CGAL {

/*!
\ingroup nt_builtin

The fundamental type `float` is an `RealEmbeddable` 
`FieldWithSqrt`. Due to rounding errors and overflow `float` is 
considered as not exact. 

\models ::FieldWithSqrt 
\models ::RealEmbeddable 

*/
class float {

}; /* end float */

/*! 
Determines whether the argument represents a value in \f$ \R\f$. 
\relates float 
*/ 
bool CGAL::is_finite(float x); 

} /* end namespace CGAL */
