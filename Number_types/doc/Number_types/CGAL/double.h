/*!
\file double.h
\ingroup nt_builtin

This header provides all necessary functions so the fundamental
type `double` is a model of the concepts `RealEmbeddable` and
`Field`. Due to rounding errors and overflow `double` is considered as
not exact.

\cgalModels `FieldWithSqrt`
\cgalModels `RealEmbeddable`
*/


namespace CGAL {
/*!
Determines whether the argument represents a value in \f$ \mathbb{R}\f$.
*/
bool is_finite(double x);

} // end namespace
