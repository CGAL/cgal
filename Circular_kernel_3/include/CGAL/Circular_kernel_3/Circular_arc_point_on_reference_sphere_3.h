#ifndef CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H
#define CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H

#include <iostream>
#include <cassert>

#include <CGAL/Circular_arc_point_3.h>

namespace CGAL {
  namespace CGALi {
  
  typedef float HQ_NT;
  
    
  template <class SK>
	class Theta_rep{
    typedef typename SK::AK::RO_t RO_t;
    typedef typename SK::FT FT;    
		typedef std::pair<HQ_NT,RO_t> Rep;
		typename CGAL_SK::template Handle<Rep>::type base;
		public:
		Theta_rep(const HQ_NT& hq,const RO_t& r):base(std::pair<HQ_NT,RO_t>(hq,r)){}
    Theta_rep(){}
    const HQ_NT& hq() const {return CGAL::get(base).first;}
    const typename SK::AK::RO_t& ftheta() const {return CGAL::get(base).second;}
	};    
    
  template<class SK>
  class Circular_arc_point_on_reference_sphere_3:public Circular_arc_point_3<SK>{
    typedef typename SK::AK::RO_t RO_t;
    typedef typename SK::FT FT;
    typedef typename CGAL_SK::Circular_arc_point_3 Circular_arc_point_3;
    //---------------
    Theta_rep<SK> Trep;
    public:
    Circular_arc_point_on_reference_sphere_3(const FT& ftheta,const FT& xt,const FT& yt,const FT& zt,const HQ_NT& _hq)
      :Circular_arc_point_3(typename CGAL_SK::Point_3(xt,yt,zt)),Trep(_hq,ftheta){};//critical point of non normal circles
        
    Circular_arc_point_on_reference_sphere_3(const HQ_NT& _hq,const RO_t& ftheta,const RO_t& x_,const RO_t& y_,const RO_t& z_)
      :Circular_arc_point_3(x_,y_,z_),Trep(_hq,ftheta){};

    Circular_arc_point_on_reference_sphere_3(const HQ_NT& _hq,const RO_t& ftheta,const typename SK::AK::Root_for_spheres_2_3& rfs)
      :Circular_arc_point_3(rfs),Trep(_hq,ftheta){};
        
    Circular_arc_point_on_reference_sphere_3():Circular_arc_point_3(){};
			
    const Theta_rep<SK>& theta_rep() const {return Trep;};
    const RO_t& get_f_of_theta() const {return theta_rep().ftheta();};  
    const HQ_NT& get_hq() const {return theta_rep().hq();}

    double get_theta_approx() const{
      double ax=CGAL::to_double(this->x());
      double ay=CGAL::to_double(this->y());
      return ( (atan2 (ay,ax)<0)?(atan2 (ay,ax)+2.*M_PI):(atan2 (ay,ax)) );
    };
    
    CGAL::Cartesian<double>::Point_3 get_point_approx() const {//just for intersection and critical points of normal circles
      return CGAL::Cartesian<double>::Point_3(CGAL::to_double(this->x()),CGAL::to_double(this->y()),CGAL::to_double(this->z()));
    }
    
    static inline Circular_arc_point_on_reference_sphere_3 VirtualPt_to_point_on_sphere(){
      return Circular_arc_point_on_reference_sphere_3(FT(0),FT(0),FT(0),FT(0),HQ_NT(-1));
    };

    static inline Circular_arc_point_on_reference_sphere_3 Root_for_sphere_to_point_on_sphere(const HQ_NT& hq,const typename AK::Root_for_spheres_2_3& R){
      return Circular_arc_point_on_reference_sphere_3(hq,auto_ftype(hq)==TAN?(R.y()/R.x()):(R.x()/R.y()),R);
    };
    
    
  };    
    
  } // namespace CGALi
} // namespace CGAL

#endif //CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H

