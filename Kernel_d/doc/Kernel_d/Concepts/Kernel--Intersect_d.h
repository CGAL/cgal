
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

A model for this must provide: 

*/

class Kernel_d::Intersect_d {
public:

/*! 
returns the result of the intersection of \f$ p\f$ and \f$ q\f$ in form of a 
polymorphic object. `Type1` and `Type2` may be any of 
`Kernel_d::Segment_d`, `Kernel_d::Ray_d`, `Kernel_d::Line_d`, 
`Kernel_d::Hyperplane_d`. 

\pre `p` and `q` have the same dimension. 
*/ 
  template <class Type1, class Type2>
boost::result_of<Kernel::Intersect_d(Type1, Type2)>::type
operator()(const Type1& p, const Type2& q); 


}; /* end Kernel_d::Intersect_d */

