#ifndef CGAL_CURVED_KERNEL_FUNCTIONS_ON_LINE_2_H
#define CGAL_CURVED_KERNEL_FUNCTIONS_ON_LINE_2_H

namespace CGAL {
namespace LinearFunctors {

  template < class CK >
  typename CK::Algebraic_1_curve_2
  get_equation( const typename CK::Line_2 & L )
  {
    typedef typename CK::RT RT;
 
    return typename CK::Algebraic_1_curve_2(RT(0),RT(0),RT(0),
					    L.a(),L.b(),L.c());
  }
  
  template < class CK >
  typename CK::Line_2  
  construct_line_2 ( const typename CK::Algebraic_1_curve_2 &eq )
  {
    return typename CK::Line_2(eq[2],eq[1],eq[0]); 
  }
  

} // namespace LinearFunctors
} // namespace CGAL
#endif // CGAL_CURVED_KERNEL_FUNCTIONS_ON_LINE_2_H
