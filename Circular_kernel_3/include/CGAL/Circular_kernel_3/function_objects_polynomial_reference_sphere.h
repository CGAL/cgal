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
#include <CGAL/Circular_kernel_3/internal_function_compare_to_left_spherical_kernel.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_reference_sphere_3.h>
#include <CGAL/Object.h>


namespace CGAL {
namespace SphericalFunctors {

  //ACCESS FUNCTIONS
  template <class SK>
  class Compute_supporting_circle_on_reference_sphere_3: Has_qrt{
    typedef typename SK::Half_circle_on_reference_sphere_3   Half_circle_on_reference_sphere_3;

  public:

    typedef typename SK::Circle_on_reference_sphere_3   result_type;
    typedef const result_type &                         qualified_result_type;
    typedef Arity_tag<1>                                Arity;

    qualified_result_type operator() (const Half_circle_on_reference_sphere_3 & a) const
    { return (a.rep().supporting_circle()); }
  };

  template <class SK>
  class Compute_half_circle_position_3: Has_qrt{
    typedef typename SK::Half_circle_on_reference_sphere_3   Half_circle_on_reference_sphere_3;

  public:

    typedef typename CGAL::Hcircle_type                 result_type;
    typedef const result_type &                         qualified_result_type;
    typedef Arity_tag<1>                                Arity;

    qualified_result_type operator() (const Half_circle_on_reference_sphere_3 & a) const
    { return (a.rep().get_position()); }
  };
  
  
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
    typedef typename SK::Circular_arc_on_reference_sphere_3 Circular_arc_on_reference_sphere_3;
  public:

    typedef typename SK::Sphere_with_radius_3 result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Circle_on_reference_sphere_3 & a) const
    { return (a.rep().reference_sphere()); }
    
    qualified_result_type operator() (const Circular_arc_on_reference_sphere_3 & a) const
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
  class Construct_circular_arc_on_reference_sphere_3
  {

    typedef typename SK::Line_3                                             Line_3;
    typedef typename SK::Point_3                                            Point_3;
    typedef typename SK::Segment_3                                          Segment_3;
    typedef typename SK::Sphere_3                                           Sphere_3;
    typedef typename SK::Plane_3                                            Plane_3;
    typedef typename SK::Line_arc_3                                         Line_arc_3;
    typedef typename SK::Circular_arc_point_on_reference_sphere_3           Circular_arc_point_on_reference_sphere_3;
    typedef typename SK::Circle_on_reference_sphere_3                       Circle_on_reference_sphere_3;
    typedef typename SK::Sphere_with_radius_3                               Sphere_with_radius_3;
    typedef typename SK::Circular_arc_on_reference_sphere_3                 Circular_arc_on_reference_sphere_3;
    typedef typename SK::Kernel_base::Circular_arc_on_reference_sphere_3    RCircular_arc_on_reference_sphere_3;
    typedef typename Circular_arc_on_reference_sphere_3::Rep                Rep;    
    
  public:
    typedef Circular_arc_on_reference_sphere_3   result_type;
    typedef Arity_tag<3> Arity; // It is not true that each constructor has
                                // 3 operands, maybe we should remove this

    result_type
    operator()(void) const
    { return Rep(); }

    result_type
    operator()(const Circle_on_reference_sphere_3 &c) const
    { return Rep(c); }

    result_type
    operator()(const Circle_on_reference_sphere_3 &l,
	       const Circular_arc_point_on_reference_sphere_3 &s,
	       const Circular_arc_point_on_reference_sphere_3 &t) const
    { return Rep(l,s,t); }

    result_type
    operator()(const Circle_on_reference_sphere_3 &l,
	       const Point_3 &s,
	       const Circular_arc_point_on_reference_sphere_3 &t) const
    { return Rep(l,s,t); }

    result_type
    operator()(const Circle_on_reference_sphere_3 &l,
	       const Circular_arc_point_on_reference_sphere_3 &s,
	       const Point_3 &t) const
    { return Rep(l,s,t); }

    result_type
    operator()(const Circle_on_reference_sphere_3 &l,
	       const Point_3 &s,
	       const Point_3 &t) const
    { return Rep(l,s,t); }

    result_type
    operator()(const Circle_on_reference_sphere_3 &c, 
               const Sphere_with_radius_3 &s1, bool less_xyz_s1,
               const Sphere_with_radius_3 &s2, bool less_xyz_s2) const
    { return Rep(c,s1,less_xyz_s1,s2,less_xyz_s2); }

    result_type
    operator()(const Sphere_with_radius_3 &s1, bool less_xyz_s1,
               const Sphere_with_radius_3 &s2, bool less_xyz_s2,
               const Circle_on_reference_sphere_3 &c) const
    { return Rep(c,s1,less_xyz_s1,s2,less_xyz_s2); }

    result_type
    operator()(const Circle_on_reference_sphere_3 &c, 
               const Plane_3 &p1, bool less_xyz_p1,
               const Plane_3 &p2, bool less_xyz_p2) const
    { return Rep(c,p1,less_xyz_p1,p2,less_xyz_p2); }

    result_type
    operator()(const Plane_3 &p1, bool less_xyz_p1,
               const Plane_3 &p2, bool less_xyz_p2,
               const Circle_on_reference_sphere_3 &c) const
    { return Rep(c,p1,less_xyz_p1,p2,less_xyz_p2); }    
  };
  
  template < class SK >
  class Construct_half_circle_on_reference_sphere_3{
    typedef typename SK::Half_circle_on_reference_sphere_3              Half_circle_on_reference_sphere_3;
    typedef typename Half_circle_on_reference_sphere_3::Rep                  Rep;
  public:
    typedef Half_circle_on_reference_sphere_3                           result_type;
    typedef Arity_tag<1>                                                Arity;

    //~ result_type
    //~ operator()()
    //~ { return Rep();}
  
    result_type
    operator()(const typename SK::Circle_on_reference_sphere_3& C,CGAL::Hcircle_type pos)
    { return Rep(C,pos);}
    
  };
  
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
    operator()()
    { return Rep();}
  
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
    
    result_type
    operator()(const CGAL::Sphere_3<SK>& S)
    { return Rep(S);}    
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

  //FUNCTORS
  template<class SK>
  class Theta_extremal_point_3{
    
    template <class OutputIterator>
    void normal_solve(const typename SK::Circle_on_reference_sphere_3& C,
                                  OutputIterator Rpts) const{
      std::vector<CGAL::Object> cont;
      typename SK::Construct_plane_3 theConstruct_plane_3 = SK().construct_plane_3_object();
      //~ typename SK::Plane_3 p=theConstruct_plane_3(typename SK::Algebraic_kernel::Polynomial_1_3(0,0,1,-C.extremal_point_z()));
      typename SK::Plane_3 p=theConstruct_plane_3(typename SK::Algebraic_kernel().construct_polynomial_1_3_object()(0,0,1,-C.extremal_point_z()));
      typename SK::Intersect_3()(C.supporting_sphere(),C.reference_sphere(),p,std::back_inserter(cont));
      CGAL_precondition(cont.size()==2);
      for (int i=0;i<2;++i)
        CGAL::assign(Rpts[i],cont[i]);
    }
    
    public:
    
    typedef void result_type;
    typedef Arity_tag< 3 >      Arity;
    
    void set_CA(const typename SK::Circle_on_reference_sphere_3& C,
                        Inter_alg_info& CA) const{
      typename std::pair<typename SK::Circular_arc_point_3,unsigned> Rpts[2];
      normal_solve(C,Rpts);
      set_IA<SK>(CA,Rpts,false);
    }
    
    template<class container>
    void set_CA(Inter_alg_info& CA,const container& Pts) const{
      set_IA<SK>(CA,Pts,false);
    }    
    
    template <class OutputIterator>
    OutputIterator operator()(const typename SK::Circle_on_reference_sphere_3& C,
                                               OutputIterator it) const{  
      Inter_alg_info CA;
      typename std::pair<typename SK::Circular_arc_point_3,unsigned> Rpts[2];
      normal_solve(C,Rpts);
      set_CA(CA,Rpts);
      *it++= typename SK::Circular_arc_point_on_reference_sphere_3(CA.qF,Rpts[CA.F_index].first);
      *it++= typename SK::Circular_arc_point_on_reference_sphere_3(CA.qS,Rpts[CA.S_index].first);
       return it;
    }
    
    template <class OutputIterator>
    OutputIterator operator()(const typename SK::Circle_on_reference_sphere_3& C,
                                               OutputIterator it,
                                               const Inter_alg_info& CA,
                                               CGAL::EvtPt_num num=CGAL::ALL) const{
      typename std::pair<typename SK::Circular_arc_point_3,unsigned> Rpts[2];
      normal_solve(C,Rpts);
      if (num!=END_PT)
        *it++= typename SK::Circular_arc_point_on_reference_sphere_3(CA.qF,Rpts[CA.F_index].first);
      else{
        *it++= typename SK::Circular_arc_point_on_reference_sphere_3(CA.qS,Rpts[CA.S_index].first);
        return it;
      }
      if (num==ALL)
        *it++= typename SK::Circular_arc_point_on_reference_sphere_3(CA.qS,Rpts[CA.S_index].first);
      return it;
    }
  };
  
  template <class SK>
  struct Compare_theta_3{
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag< 2 >          Arity;
    
    result_type operator()(const typename SK::Theta_rep& T1,const typename SK::Theta_rep& T2) const {
      CGAL::HQ_NT res=T1.hq()-T2.hq();
      if (res == 0.){//same quadrant 
        int m=(CGAL::auto_ftype(T1.hq())==CGAL::TAN)?(1):(-1);
        if (truncf(T1.hq())!=T1.hq())//for hquadrant boundary
          res=0;
        else
          res= m*CGAL::compare(T1.ftheta(),T2.ftheta());
      }
      if (res < 0) return (CGAL::SMALLER);
      if (res >0 ) return (CGAL::LARGER);  
      return CGAL::EQUAL;
    }

    result_type operator()(const typename SK::Circular_arc_point_on_reference_sphere_3& p1,
                                       const typename SK::Circular_arc_point_on_reference_sphere_3& p2) const{
      return (*this)(p1.theta_rep(),p2.theta_rep());
    }    
    
  };

  template <class SK>
  struct Compare_theta_z_3{
    typedef CGAL::Comparison_result   result_type;
    typedef Arity_tag< 2 >            Arity;
    
    result_type operator()(const typename CGAL::Circular_arc_point_on_reference_sphere_3<SK>& p1,
                           const typename CGAL::Circular_arc_point_on_reference_sphere_3<SK>& p2,
                           bool decreasing_z=false) const{
      SK sk;
      CGAL::Comparison_result res=sk.compare_theta_3_object()(p1,p2);
      if (res==CGAL::EQUAL)
        res=decreasing_z?(CGAL::opposite(sk.compare_z_3_object()(p1,p2))):(sk.compare_z_3_object()(p1,p2));
      return res;
    }
  };
  
  template <class SK>
  class Compare_z_at_theta_3{
    
    typedef typename SK::Half_circle_on_reference_sphere_3 Half_circle_on_reference_sphere_3;
    typedef typename SK::Circle_on_reference_sphere_3 Circle_on_reference_sphere_3;
    typedef typename SK::Circular_arc_point_on_reference_sphere_3 Circular_arc_point_on_reference_sphere_3;
    typedef typename SK::FT FT;
    
    public:
    typedef CGAL::Comparison_result   result_type;
    typedef Arity_tag< 2 >            Arity;
    
    //At theta=0
    //Compare the intersection point of an half circle with M_theta for theta=0 with the critical points of a circle on reference sphere
    result_type operator() (const Half_circle_on_reference_sphere_3& H,
                                            const Circle_on_reference_sphere_3& S) const{
      std::vector<CGAL::Object> cont;
      SK sk;
      typename SK::Plane_3 p=sk.construct_plane_3_object()(typename SK::Algebraic_kernel::Polynomial_1_3(0,1,0,0));
      typename SK::Circular_arc_point_3 Pt;
      sk.intersect_3_object()(H.supporting_circle().reference_sphere(),H.supporting_circle().supporting_sphere(),p,std::back_inserter(cont));
      select_inter_pt(H.supporting_circle().type_of_circle_on_reference_sphere()==CGAL::THREADED?CGAL::UNDEF:H.get_position(),cont,Pt,0.5);      
      return sk.compare_z_3_object()(Pt,typename SK::Point_3(0,0,S.extremal_point_z()));
    }
    //theta=0
    
    //Compare the z coordinates of the intersection points with the meridian at theta=0 of two half_circles
    result_type operator() (const Half_circle_on_reference_sphere_3 & H1,
                                            const Half_circle_on_reference_sphere_3& H2) const{
      return (*this)(H1,H2,0.5,0);
    }
    
    //Compare the z coordinates of the intersection points with the meridian defined by hq and A (A=f(theta) compliant with hq) of two half circles.
    result_type operator() ( const Half_circle_on_reference_sphere_3& H1,
                                            const Half_circle_on_reference_sphere_3 & H2,
                                            const typename CGAL::HQ_NT& hq, const FT& A) const {
      std::vector<CGAL::Object> cont;
      FT a,b;
      if (truncf(hq)!=hq){//for hquadrant boundary
        if (hq==1.5 || hq==5.5){a=1;b=-1;}
        else if (hq==3.5 || hq==7.5){a=1;b=1;}
        else if (hq==0.5 || hq==4.5 || hq==8.5){a=0;b=1;}
        else if (hq==2.5 || hq==6.5){a=1;b=0;}
      }
      else
        if (CGAL::auto_ftype(hq)==CGAL::TAN){a=A;b=-1;}
        else{a=-1;b=A;}
      
      SK sk;
      typename SK::Plane_3 p=sk.construct_plane_3_object()(typename SK::Algebraic_kernel::Polynomial_1_3(a,b,0,0));
      typename SK::Circular_arc_point_3 Pt1,Pt2;
      sk.intersect_3_object()(H1.supporting_circle().reference_sphere(),H1.supporting_circle().supporting_sphere(),p,std::back_inserter(cont));
      select_inter_pt(H1.supporting_circle().type_of_circle_on_reference_sphere()==CGAL::THREADED?CGAL::UNDEF:H1.get_position(),cont,Pt1,hq);
      cont.clear();
      sk.intersect_3_object()(H2.supporting_circle().reference_sphere(),H2.supporting_circle().supporting_sphere(),p,std::back_inserter(cont));
      select_inter_pt(H2.supporting_circle().type_of_circle_on_reference_sphere()==CGAL::THREADED?CGAL::UNDEF:H2.get_position(),cont,Pt2,hq);
      return sk.compare_z_3_object()(Pt1,Pt2);
    }

    //compare a normal start point vs a theta-monotonic circle arc.
    result_type operator() (const Circular_arc_point_on_reference_sphere_3& Pt,
                                            const Half_circle_on_reference_sphere_3& H) const{
      CGAL::Hcircle_type pos=H.get_position();
      if ( (pos==CGAL::SENT_SPOLE) || (pos==CGAL::SENT_NPOLE) )
      {
        return (pos==CGAL::SENT_SPOLE)?(CGAL::LARGER):(CGAL::SMALLER);
      }
      if (H.supporting_circle().type_of_circle_on_reference_sphere()==CGAL::NORMAL)
      {
        int center_pos=0;//will contain Pt.z-start_pt(H).z
        int res = point_VS_supporting_plane(Pt,H.supporting_circle(),center_pos);
        if (res > 0)
          return (pos==CGAL::UPPER)?(CGAL::SMALLER):(CGAL::LARGER);
        if (res < 0 )
          return (center_pos>0)?(CGAL::LARGER):(CGAL::SMALLER);
        //case res=0
        if (pos==CGAL::UPPER)
          return (center_pos>0)?(CGAL::EQUAL):(CGAL::SMALLER);
        return (center_pos<0)?(CGAL::EQUAL):(CGAL::LARGER);
      }
      int center_pos=-1;
      int res = point_VS_supporting_plane(Pt,H.supporting_circle(),center_pos);
      if (res==0)
        return (CGAL::EQUAL);
      res*=signof(H.supporting_circle().supporting_sphere_center().z());//normal vector always to the top
      return (res<0)?(CGAL::SMALLER):(CGAL::LARGER);//threaded or polar circle
    }

    private:
    //Select the intersection point of meridian with Half  circle
    void select_inter_pt(CGAL::Hcircle_type T,
                                    const std::vector<CGAL::Object>& cont,
                                    typename SK::Circular_arc_point_3& Pt,
                                    const typename CGAL::HQ_NT& hq) const{
      std::pair<typename SK::Circular_arc_point_3,unsigned> p1,p2;
      CGAL::assign(p1,cont[0]);
      CGAL::assign(p2,cont[1]);
      SK sk;
      if (T==CGAL::UNDEF){
        if (hq<2 || hq>7){
          Pt=(sk.compare_x_3_object()(p1.first,p2.first)==CGAL::SMALLER?p2.first:p1.first);
          return;
        }
        if (hq>=2 && hq <=3){
          Pt=(sk.compare_y_3_object()(p1.first,p2.first)==CGAL::SMALLER?p2.first:p1.first);
          return;
        }
        if (hq>3 && hq<6){
          Pt=(sk.compare_x_3_object()(p1.first,p2.first)==CGAL::SMALLER?p1.first:p2.first );
          return;
        }
        if (hq>=6 || hq <=7){
          Pt=(sk.compare_y_3_object()(p1.first,p2.first)==CGAL::SMALLER?p1.first:p2.first);
          return;
        }        
      }
      else{
        if (sk.compare_z_3_object()(p1.first,p2.first)==CGAL::SMALLER)
          Pt=(T==CGAL::LOWER?p1.first:p2.first);
        else
          Pt=(T==CGAL::LOWER?p2.first:p1.first);
      }
    }

    //return 1,0,-1  if up to ,on ,under to the plane (else for threaded circle X sign_of(z) to have the same result)
     int point_VS_supporting_plane(const Circular_arc_point_on_reference_sphere_3& pt,
                                                     const Circle_on_reference_sphere_3& C,
                                                     int& center_pos) const{
      typename SK::Plane_3 plane=C.supporting_plane();
      int i=typename SK::Algebraic_kernel().sign_at_object()(SK().get_equation_object()(plane),pt.coordinates());
      if (center_pos==0)//compute the position of the SP wrt circle center only for normal circles
      {
        center_pos=CGAL::sign(pt.z()-C.extremal_point_z());
        i*=signof(C.circle_center_coefficient());//signof(...) handle IVM
      }
      return i;
    }
  };  
  
  template <class SK>
  struct Compare_z_to_left_3{
    typedef CGAL::Comparison_result      result_type;
    typedef Arity_tag< 2 >                      Arity;
    //theta=0
    result_type 
    operator()( const typename SK::Half_circle_on_reference_sphere_3& H1,
                          const typename SK::Half_circle_on_reference_sphere_3& H2) const {
      return CGAL::SphericalFunctors::compare_to_hcircle_to_left<SK,CGAL::SphericalFunctors::trait_for_cmp_tgt_theta_0<SK> >()(H1,H2);
    }
    
    result_type 
    operator()( const typename SK::Half_circle_on_reference_sphere_3& H1,
                          const typename SK::Half_circle_on_reference_sphere_3& H2,
                          const typename SK::Circular_arc_point_on_reference_sphere_3& p) const{
      CGAL::SphericalFunctors::compare_to_hcircle_to_left<SK,CGAL::SphericalFunctors::trait_for_cmp_tgt<SK> > cmp(p.coordinates());
      return cmp(H1,H2);
    }
  };  
  
  
}
}

#endif
