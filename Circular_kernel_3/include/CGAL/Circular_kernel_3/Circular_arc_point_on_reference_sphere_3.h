#ifndef CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H
#define CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H

#include <iostream>
#include <cassert>

#include <CGAL/Circular_arc_point_3.h>


//~ enum Fct_type{TAN, COT, FIXED, TAG_M2};

//~ inline Fct_type auto_ftype(const float& hquad){//SHOULD BE HQ_NT
//~ if (hquad >7 || hquad<2 || (hquad>3 && hquad<6))
  //~ return TAN;
//~ return COT;
//~ };


namespace CGAL {
  namespace CGALi {
  
    
  template<class SK>
  class Circular_arc_point_on_reference_sphere_3:public Circular_arc_point_3<SK>{
    typedef typename SK::Algebraic_kernel::Root_of_2 Root_of_2;
    typedef typename SK::FT FT;
    typedef typename SK::Algebraic_kernel AK;
    typedef Circular_arc_point_3<SK> T_Circular_arc_point_3;
    typedef typename SK::HQ_NT HQ_NT;
    typedef typename SK::Theta_rep Theta_rep;
    //---------------
    Theta_rep Trep;
    public:
    Circular_arc_point_on_reference_sphere_3(const FT& ftheta,const FT& xt,const FT& yt,const FT& zt,const HQ_NT& _hq)
      :T_Circular_arc_point_3(typename SK::Point_3(xt,yt,zt)),Trep(_hq,ftheta){};//critical point of non normal circles
        
    Circular_arc_point_on_reference_sphere_3(const HQ_NT& _hq,const Root_of_2& ftheta,const Root_of_2& x_,const Root_of_2& y_,const Root_of_2& z_)
      :T_Circular_arc_point_3(x_,y_,z_),Trep(_hq,ftheta){};

    Circular_arc_point_on_reference_sphere_3(const HQ_NT& _hq,const Root_of_2& ftheta,const typename SK::Algebraic_kernel::Root_for_spheres_2_3& rfs)
      :T_Circular_arc_point_3(rfs),Trep(_hq,ftheta){};
        
    Circular_arc_point_on_reference_sphere_3():T_Circular_arc_point_3(FT(0),FT(0),FT(0)),Trep(-1,FT(0)){};
      
    Circular_arc_point_on_reference_sphere_3(const HQ_NT& hq,const typename AK::Root_for_spheres_2_3& R):T_Circular_arc_point_3(R),Trep(hq,auto_ftype(hq)==TAN?(R.y()/R.x()):(R.x()/R.y())){};            
      
    //~ static inline Circular_arc_point_on_reference_sphere_3 VirtualPt_to_point_on_sphere(){
      //~ return Circular_arc_point_on_reference_sphere_3(FT(0),FT(0),FT(0),FT(0),HQ_NT(-1));
    //~ };

    //~ static inline Circular_arc_point_on_reference_sphere_3 Root_for_sphere_to_point_on_sphere(const HQ_NT& hq,const typename AK::Root_for_spheres_2_3& R){
      //~ return Circular_arc_point_on_reference_sphere_3(hq,auto_ftype(hq)==TAN?(R.y()/R.x()):(R.x()/R.y()),R);
    //~ };      
      
      
    const Theta_rep& theta_rep() const {return Trep;};
    const Root_of_2& get_f_of_theta() const {return theta_rep().ftheta();};  
    const HQ_NT& get_hq() const {return theta_rep().hq();}

    double get_theta_approx() const{
      double ax=CGAL::to_double(this->x());
      double ay=CGAL::to_double(this->y());
      return ( (atan2 (ay,ax)<0)?(atan2 (ay,ax)+2.*M_PI):(atan2 (ay,ax)) );
    };
    
    CGAL::Cartesian<double>::Point_3 get_point_approx() const {//just for intersection and critical points of normal circles
      return CGAL::Cartesian<double>::Point_3(CGAL::to_double(this->x()),CGAL::to_double(this->y()),CGAL::to_double(this->z()));
    }
    
    
    
  };    
    
  } // namespace CGALi
} // namespace CGAL

#endif //CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H

