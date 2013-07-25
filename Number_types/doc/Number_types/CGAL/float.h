

/*!
\file float.h
\ingroup nt_builtin

This header provides all necessary functions so the fundamental type
`float` is a model of the concepts `RealEmbeddable` and
`FieldWithSqrt`. Due to rounding errors and overflow `float` is
considered as not exact.

\cgalModels `FieldWithSqrt` 
\cgalModels `RealEmbeddable` 

*/

namespace CGAL {

/*!
Determines whether the argument represents a value in \f$ \mathbb{R}\f$. 
*/ 
bool is_finite(float x); 

} /* end namespace CGAL */
