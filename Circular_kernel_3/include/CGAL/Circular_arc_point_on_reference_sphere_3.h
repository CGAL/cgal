#ifndef CGAL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H
#define CGAL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H

namespace CGAL {

template < typename SphericalKernel >
class Circular_arc_point_on_reference_sphere_3
  : public SphericalKernel::Kernel_base::Circular_arc_point_on_reference_sphere_3
{
  typedef typename SphericalKernel::Kernel_base::Circular_arc_point_on_reference_sphere_3 
                                           RCircular_arc_point_on_reference_sphere_3;

  typedef typename SphericalKernel::Root_of_2             Root_of_2;
  typedef typename SphericalKernel::Point_3               Point_3;
  typedef typename SphericalKernel::Plane_3               Plane_3;
  typedef typename SphericalKernel::Line_3                Line_3;
  typedef typename SphericalKernel::Circle_3              Circle_3;
  typedef typename SphericalKernel::Sphere_3              Sphere_3;
  typedef typename SphericalKernel::FT                    FT;
  typedef typename SphericalKernel::Algebraic_kernel      AK;
  typedef typename SphericalKernel::HQ_NT                 HQ_NT;

public:
  typedef typename SphericalKernel::Root_for_spheres_2_3 
    Root_for_spheres_2_3;
  typedef SphericalKernel   R; 
  typedef RCircular_arc_point_on_reference_sphere_3 Repd;

  const Repd& rep() const
  {
    return *this;
  }

  Repd& rep()
  {
    return *this;
  }

  Circular_arc_point_on_reference_sphere_3(const RCircular_arc_point_on_reference_sphere_3& p)
  :RCircular_arc_point_on_reference_sphere_3(p){}    
  
  Circular_arc_point_on_reference_sphere_3()
  : RCircular_arc_point_on_reference_sphere_3(
    typename R::Construct_circular_arc_point_on_reference_sphere_3()())
  {}

  Circular_arc_point_on_reference_sphere_3(const FT& ftheta,const FT& xt,const FT& yt,const FT& zt,const HQ_NT& _hq)
  : RCircular_arc_point_on_reference_sphere_3(
    typename R::Construct_circular_arc_point_on_reference_sphere_3()(ftheta,xt,yt,zt,_hq))
  {}
    
  Circular_arc_point_on_reference_sphere_3(const HQ_NT& _hq,const Root_of_2& ftheta,const Root_of_2& x_,const Root_of_2& y_,const Root_of_2& z_)
  : RCircular_arc_point_on_reference_sphere_3(
    typename R::Construct_circular_arc_point_on_reference_sphere_3()(_hq,ftheta,x_,y_,z_))
  {}
    
  Circular_arc_point_on_reference_sphere_3(const HQ_NT& hq,const typename AK::Root_for_spheres_2_3& R)
  : RCircular_arc_point_on_reference_sphere_3(hq,R)
  {}
    
  Circular_arc_point_on_reference_sphere_3(const HQ_NT& hq,const typename SphericalKernel::Circular_arc_point_3& R)
  : RCircular_arc_point_on_reference_sphere_3(hq,R)
  {}    

  Circular_arc_point_on_reference_sphere_3(const HQ_NT& _hq,const Root_of_2& ftheta,const typename AK::Root_for_spheres_2_3& rfs)
  : RCircular_arc_point_on_reference_sphere_3(
    typename R::Construct_circular_arc_point_on_reference_sphere_3()(_hq,ftheta,rfs))
  {}
    

  typename Qualified_result_of<typename R::Compute_circular_theta_rep_3,Circular_arc_point_on_reference_sphere_3>::type
  theta_rep() const
  { return typename R::Compute_circular_theta_rep_3()(*this);}
  
  
      
  typename Qualified_result_of<typename R::Compute_circular_theta_3,Circular_arc_point_on_reference_sphere_3>::type
  //const Root_of_2 &
  get_f_of_theta() const
  { return typename R::Compute_circular_theta_3()(*this);}  

  typename Qualified_result_of<typename R::Compute_circular_hq_3,Circular_arc_point_on_reference_sphere_3>::type
  get_hq() const
  { return typename R::Compute_circular_hq_3()(*this);}    
  
  //~ Bbox_3 bbox() const
  //~ { return typename R::Construct_bbox_3()(*this); }

};


}

#endif
