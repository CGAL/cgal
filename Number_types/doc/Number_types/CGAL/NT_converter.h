namespace CGAL {

/*!
\ingroup nt_util

A number type converter usable as default, for CGAL::Cartesian_converter and CGAL::Homogeneous_converter.

\models ::AdaptableFunctor 

*/
template < class NT1, class NT2 >
class NT_converter{
public:

/// \name Operations
///@{
/*! 
convert `a` from NT1 to NT2.
*/ 
NT2 operator()(const NT1 &a) const;

/// @}

}; /* end NT_converter */
} /* end namespace CGAL */