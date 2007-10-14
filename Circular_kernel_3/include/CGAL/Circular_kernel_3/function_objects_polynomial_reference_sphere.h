#ifndef CGAL_SPHERICAL_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_REFERENCE_SPHERE_H
#define CGAL_SPHERICAL_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_REFERENCE_SPHERE_H

#include <CGAL/kernel_basic.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_circular_arc_point_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_sphere_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_line_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_plane_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_circle_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_line_arc_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_circular_arc_3.h>
#include <CGAL/Circular_kernel_3/internal_function_has_on_spherical_kernel.h>
#include <CGAL/Circular_kernel_3/internal_function_compare_spherical_kernel.h>
#include <CGAL/Object.h>


namespace CGAL {
namespace SphericalFunctors {

  //ACCESS FUNCTIONS
  template <class SK>
  class Compute_circle_center_coefficient_3/*: Has_qrt*/{
    typedef typename SK::Circle_on_reference_sphere_3   Circle_on_reference_sphere_3;

  public:

    typedef typename SK::FT result_type;
    //~ typedef const result_type &        qualified_result_type;
    typedef result_type               qualified_result_type;
    typedef Arity_tag<1>              Arity;

    qualified_result_type operator() (const Circle_on_reference_sphere_3 & a) const
    { return (a.rep().circle_center_coefficient()); }
  };
  
  template <class SK>
  class Compute_extremal_point_z/*: Has_qrt*/{
    typedef typename SK::Circle_on_reference_sphere_3   Circle_on_reference_sphere_3;

  public:

    typedef typename SK::FT result_type;
    //~ typedef const result_type &        qualified_result_type;
    typedef result_type               qualified_result_type;
    typedef Arity_tag<1>              Arity;

    qualified_result_type operator() (const Circle_on_reference_sphere_3 & a) const
    { return (a.rep().extremal_point_z()); }
  };
  
  
  template <class SK>
  class Compute_type_of_circle_on_reference_sphere_3: Has_qrt{
    typedef typename SK::Circle_on_reference_sphere_3   Circle_on_reference_sphere_3;

  public:

    typedef typename CGAL::Circle_type result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Circle_on_reference_sphere_3 & a) const
    { return (a.rep().type_of_circle_on_reference_sphere()); }    
  };
  
  template <class SK>
  class Compute_supporting_sphere_radius_3: Has_qrt{
    typedef typename SK::Circle_on_reference_sphere_3   Circle_on_reference_sphere_3;

  public:

    typedef typename SK::FT result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Circle_on_reference_sphere_3 & a) const
    { return (a.rep().supporting_sphere_radius()); }    
  };
  
  //BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD 
  template <class SK>
  class Compute_supporting_sphere_squared_radius_3/*: Has_qrt*/{
    typedef typename SK::Circle_on_reference_sphere_3   Circle_on_reference_sphere_3;

  public:

    typedef typename SK::FT result_type;
    //~ typedef const result_type &        qualified_result_type;
    typedef result_type        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Circle_on_reference_sphere_3 & a) const
    { return (a.rep().supporting_sphere_squared_radius()); }    
  };  
  
  //BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD 
  template <class SK>
  class Compute_supporting_sphere_center_3/*: Has_qrt*/{
    typedef typename SK::Circle_on_reference_sphere_3   Circle_on_reference_sphere_3;

  public:

    typedef typename SK::Point_3 result_type;
    //~ typedef const result_type &        qualified_result_type;
    typedef result_type        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Circle_on_reference_sphere_3 & a) const
    { return (a.rep().supporting_sphere_center()); }
  };
  
  template <class SK>
  class Compute_reference_sphere_3: Has_qrt{
    typedef typename SK::Circle_on_reference_sphere_3   Circle_on_reference_sphere_3;

  public:

    typedef typename SK::Sphere_with_radius_3 result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Circle_on_reference_sphere_3 & a) const
    { return (a.rep().reference_sphere()); }    
  };  
  
  
  
  
  template <class SK>
  class Compute_radius_sphere_with_radius_3: Has_qrt{
    typedef typename SK::Sphere_with_radius_3   Sphere_with_radius_3;

  public:

    typedef typename SK::FT result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Sphere_with_radius_3 & a) const
    { return (a.rep().radius()); }    
  };
    
  template <class SK>
  class Compute_theta_hq_3: Has_qrt{
    typedef typename SK::Theta_rep   Theta_rep;

  public:

    typedef CGAL::HQ_NT result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Theta_rep & a) const
    { return (a.rep().hq()); }    
  };

  template <class SK>
  class Compute_theta_ftheta_3: Has_qrt{
    typedef typename SK::Theta_rep   Theta_rep;
    typedef typename SK::Root_of_2  Root_of_2;
  public:

    typedef  Root_of_2  result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Theta_rep & a) const
    { return (a.rep().ftheta()); }    
  };
  
  
  template <class SK>
  class Compute_circular_theta_rep_3: Has_qrt{
    typedef typename SK::Circular_arc_point_on_reference_sphere_3   Circular_arc_point_on_reference_sphere_3;

  public:

    typedef  typename  SK::Theta_rep result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Circular_arc_point_on_reference_sphere_3 & a) const
    { return (a.rep().theta_rep()); }
  };
  
  template <class SK>
  class Compute_circular_theta_3: Has_qrt{
    typedef typename SK::Circular_arc_point_on_reference_sphere_3   Circular_arc_point_on_reference_sphere_3;
    typedef typename SK::Root_of_2                 Root_of_2;

  public:

    typedef  Root_of_2  result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Circular_arc_point_on_reference_sphere_3 & a) const
    { return (a.rep().theta_rep().ftheta()); }    
  };

  template <class SK>
  class Compute_circular_hq_3: Has_qrt{
    typedef typename SK::Circular_arc_point_on_reference_sphere_3   Circular_arc_point_on_reference_sphere_3;

  public:

    typedef CGAL::HQ_NT result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Circular_arc_point_on_reference_sphere_3 & a) const
    { return (a.rep().theta_rep().hq()); }    
  };  
  
  //CONSTRUCTIONS

  template < class SK >
  class Construct_circle_on_reference_sphere_3
  {
    typedef typename SK::Circle_on_reference_sphere_3                           Circle_on_reference_sphere_3;
    typedef typename SK::Kernel_base::Circle_on_reference_sphere_3              RCircle_on_reference_sphere_3;
    typedef typename Circle_on_reference_sphere_3::Repd                         Rep;

  public:
    typedef  Circle_on_reference_sphere_3 result_type;
    typedef Arity_tag<1>             Arity;


    result_type
    operator()(const typename SK::FT& _r,const typename SK::Point_3& _c,const typename SK::Sphere_with_radius_3& S)
    { return Rep(_r,_c,S);}
    
    result_type
    operator()(const typename SK::FT& _r,const typename SK::Point_3& _c,CGAL::Circle_type nat,const typename SK::Sphere_with_radius_3& S)
    { return Rep(_r,_c,nat,S);}    
  };
    
  template < class SK >
  class Construct_sphere_with_radius_3
  {
    typedef typename SK::Sphere_with_radius_3                           Sphere_with_radius_3;
    typedef typename SK::Kernel_base::Sphere_with_radius_3              RSphere_with_radius_3;
    typedef typename Sphere_with_radius_3::Repd                         Rep;

  public:
    typedef  Sphere_with_radius_3 result_type;
    typedef Arity_tag<1>             Arity;


    result_type
    operator()(void) 
    { return Rep(); }
    
    result_type
    operator()(const typename SK::FT& _r,const typename SK::Point_3& _c)
    { return Rep(_r,_c);}
  };

      
  template < class SK >
  class Construct_theta_rep
  {
    typedef typename SK::Theta_rep                           Theta_rep;
    typedef typename SK::Kernel_base::Theta_rep      RTheta_rep;
    typedef typename SK::Root_of_2                            Root_of_2;
    typedef typename Theta_rep::Rep                         Rep;


  public:
    typedef  Theta_rep result_type;
    typedef Arity_tag<1>             Arity;


    result_type
    operator()(void) 
    { return Rep(); }
    
    result_type
    operator()(const HQ_NT& hq,const Root_of_2& r)
    { return Rep(hq,r);}
  };    
  
  template < class SK >
  class Construct_circular_arc_point_on_reference_sphere_3
  {
    typedef typename SK::Point_3                            Point_3;
    typedef typename SK::Plane_3                            Plane_3;
    typedef typename SK::Line_3                             Line_3;
    typedef typename SK::Circle_3                           Circle_3;
    typedef typename SK::Sphere_3                           Sphere_3;
    typedef typename SK::Circular_arc_point_on_reference_sphere_3               Circular_arc_point_on_reference_sphere_3;
    typedef typename SK::Kernel_base::Circular_arc_point_on_reference_sphere_3  RCircular_arc_point_on_reference_sphere_3;
    typedef typename SK::Root_of_2                          Root_of_2;
    typedef typename SK::FT                                     FT;
    typedef typename Circular_arc_point_on_reference_sphere_3::Repd              Rep;



  public:
    typedef  Circular_arc_point_on_reference_sphere_3 result_type;
    typedef Arity_tag<1>             Arity;


    result_type
    operator()(void) 
    { return Rep(); }
    
   
    result_type
    operator()(const FT& ftheta, const FT& xt, const FT& yt, const FT& zt, const HQ_NT& hq)
    { return Rep(ftheta,xt,yt,zt,hq);}
    
        
    result_type
    operator()(const HQ_NT& _hq,const Root_of_2& ftheta,const Root_of_2& x_,const Root_of_2& y_,const Root_of_2& z_)
    {Rep(_hq,ftheta,x_,y_,z_);}

    result_type
    operator()(const HQ_NT& _hq,const Root_of_2& ftheta,const typename SK::Algebraic_kernel::Root_for_spheres_2_3& rfs)
    {return Rep(_hq,ftheta,rfs);}

    result_type
    operator()(const HQ_NT& hq,const typename SK::Algebraic_kernel::Root_for_spheres_2_3& R)
    {return Rep(hq,R);}

    result_type
    operator()(const HQ_NT& hq,const typename SK::Circular_arc_point_3& R)
    {return Rep(hq,R);}    
  };  
  
}
}

#endif
