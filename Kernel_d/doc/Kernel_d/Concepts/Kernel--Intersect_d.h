
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Intersect_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the result of the intersection of \f$ p\f$ and \f$ q\f$ in form of a 
stack-based discriminated union container object. `Type1` and `Type2` may be any of 
`Kernel_d::Segment_d`, `Kernel_d::Ray_d`, `Kernel_d::Line_d`, 
`Kernel_d::Hyperplane_d`. 

For a list of the possible return types, see `CGAL::intersection()`.

\pre `p` and `q` have the same dimension. 
*/ 
template <class Type1, class Type2>
cpp11::result_of<Kernel::Intersect_d(Type1, Type2)>::type
operator()(const Type1& p, const Type2& q); 

/// @}

}; /* end Kernel_d::Intersect_d */

