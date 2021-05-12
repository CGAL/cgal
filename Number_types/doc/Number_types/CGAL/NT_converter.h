namespace CGAL {

/*!
\ingroup nt_util

A number type converter usable as default, for `Cartesian_converter` and `Homogeneous_converter`.

\cgalModels `AdaptableFunctor`

*/
template < class NT1, class NT2 >
struct NT_converter{

/// \name Operations
///@{
/*!
convert `a` from NT1 to NT2.
*/
NT2 operator()(const NT1 &a) const;

/// @}

}; /* end NT_converter */
} /* end namespace CGAL */
