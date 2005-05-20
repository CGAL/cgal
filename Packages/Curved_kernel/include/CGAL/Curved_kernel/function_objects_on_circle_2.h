#ifndef CGAL_CURVED_KERNEL_FUNCTION_OBJECTS_ON_CIRCLE_2_H
#define CGAL_CURVED_KERNEL_FUNCTION_OBJECTS_ON_CIRCLE_2_H

#include <CGAL/Curved_kernel/internal_functions_on_circle_2.h>

#include <CGAL/Curved_kernel/function_objects_on_line_2.h> 
// to be removed when CGAL::Kernel has a Get_equation

namespace CGAL {
namespace CircularFunctors {

  template < class CK >
  class Construct_circle_2 : public  CK::Linear_kernel::Construct_circle_2
  {
    public:

    typedef typename CK::Circle_2 result_type;

    result_type
    operator() ( const typename CK::Polynomial_for_circles_2_2 &eq )
      {
	return construct_circle_2<CK>(eq);
      }
  };

  template < class CK >
  class Get_equation : public LinearFunctors::Get_equation<CK>
  {
    public:

    typedef typename CK::Polynomial_for_circles_2_2 result_type;

    using LinearFunctors::Get_equation<CK>::operator();

    result_type
    operator() ( const typename CK::Circle_2 & c )
      {
	return get_equation<CK>(c);
      }
  };

} // namespace CircularFunctors
} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_FUNCTION_OBJECTS_ON_CIRCLE_2_H
