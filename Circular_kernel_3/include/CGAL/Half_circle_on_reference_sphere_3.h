#ifndef CGAL_HALF_CIRCLE_ON_REFERENCE_SPHERE_3_H
#define CGAL_HALF_CIRCLE_ON_REFERENCE_SPHERE_3_H

namespace CGAL {
  
template < typename SphericalKernel >
class Half_circle_on_reference_sphere_3
  : public SphericalKernel::Kernel_base::Half_circle_on_reference_sphere_3
{
  typedef typename SphericalKernel::Kernel_base::Half_circle_on_reference_sphere_3 
                                           RHalf_circle_on_reference_sphere_3;

public:
  typedef SphericalKernel   R; 
  typedef RHalf_circle_on_reference_sphere_3 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  //~ Half_circle_on_reference_sphere_3()
  //~ :RHalf_circle_on_reference_sphere_3(){}
  
  Half_circle_on_reference_sphere_3(const RHalf_circle_on_reference_sphere_3& p)
  :RHalf_circle_on_reference_sphere_3(p){}

  Half_circle_on_reference_sphere_3(const typename SphericalKernel::Circle_on_reference_sphere_3& C,CGAL::Hcircle_type pos)    
  : RHalf_circle_on_reference_sphere_3(typename R::Construct_half_circle_on_reference_sphere_3()(C,pos))
  {}

  typename Qualified_result_of<typename R::Compute_half_circle_position_3,Half_circle_on_reference_sphere_3>::type
  get_position() const
  { return typename R::Compute_half_circle_position_3()(*this);}
    
  typename Qualified_result_of<typename R::Compute_supporting_circle_on_reference_sphere_3,Half_circle_on_reference_sphere_3>::type
  supporting_circle() const
  { return typename R::Compute_supporting_circle_on_reference_sphere_3()(*this);}
   
    
};
  
}
#endif

