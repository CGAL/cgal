#ifndef CGAL_CIRCLE_ON_REFERENCE_SPHERE_3_H
#define CGAL_CIRCLE_ON_REFERENCE_SPHERE_3_H

namespace CGAL {

template < typename SphericalKernel >
class Circle_on_reference_sphere_3
  : public SphericalKernel::Kernel_base::Circle_on_reference_sphere_3
{
  typedef typename SphericalKernel::Kernel_base::Circle_on_reference_sphere_3 
                                           RCircle_on_reference_sphere_3;

  typedef typename SphericalKernel::Point_3               Point_3;
  typedef typename SphericalKernel::FT                    FT;
  typedef typename SphericalKernel::Algebraic_kernel      AK;
  typedef typename SphericalKernel::Root_for_spheres_2_3  Root_for_spheres_2_3;  

public:
  typedef SphericalKernel   R; 
  typedef RCircle_on_reference_sphere_3 Repd;

  const Repd& rep() const
  {
    return *this;
  }

  Repd& rep()
  {
    return *this;
  }

  Circle_on_reference_sphere_3(const RCircle_on_reference_sphere_3& p)
  :RCircle_on_reference_sphere_3(p){}

  Circle_on_reference_sphere_3(const FT& _r,const Point_3& _c,const typename SphericalKernel::Sphere_with_radius_3& S)
  : RCircle_on_reference_sphere_3(
    typename R::Construct_circle_on_reference_sphere_3()(_r,_c,S))
  {}

  Circle_on_reference_sphere_3(const FT& _r,const Point_3& _c,CGAL::Circle_type nat,const typename SphericalKernel::Sphere_with_radius_3& S)
  : RCircle_on_reference_sphere_3(
    typename R::Construct_circle_on_reference_sphere_3()(_r,_c,nat,S))
  {}    


  typename Qualified_result_of<typename R::Compute_type_of_circle_on_reference_sphere_3,Circle_on_reference_sphere_3>::type
  type_of_circle_on_reference_sphere() const
  { return typename R::Compute_type_of_circle_on_reference_sphere_3()(*this);}
    
  typename Qualified_result_of<typename R::Compute_supporting_sphere_radius_3,Circle_on_reference_sphere_3>::type
  supporting_sphere_radius() const
  { return typename R::Compute_supporting_sphere_radius_3()(*this);}
  
  typename Qualified_result_of<typename R::Compute_supporting_sphere_squared_radius_3,Circle_on_reference_sphere_3>::type
  supporting_sphere_squared_radius() const
  { return typename R::Compute_supporting_sphere_squared_radius_3()(*this);}
    
  typename Qualified_result_of<typename R::Compute_supporting_sphere_center_3,Circle_on_reference_sphere_3>::type
  supporting_sphere_center() const
  { return typename R::Compute_supporting_sphere_center_3()(*this);}
  
  typename Qualified_result_of<typename R::Compute_reference_sphere_3,Circle_on_reference_sphere_3>::type
  reference_sphere() const
  { return typename R::Compute_reference_sphere_3()(*this);}
  
  typename Qualified_result_of<typename R::Compute_extremal_point_z,Circle_on_reference_sphere_3>::type
  extremal_point_z() const
  { return typename R::Compute_extremal_point_z()(*this);} 
  
  typename Qualified_result_of<typename R::Compute_circle_center_coefficient_3,Circle_on_reference_sphere_3>::type
  circle_center_coefficient() const
  { return typename R::Compute_circle_center_coefficient_3()(*this);}  
  
   
    
};

}
#endif

