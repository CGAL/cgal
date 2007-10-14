#ifndef CGAL_SPHERICAL_KERNEL_CONSTANT_H
#define CGAL_SPHERICAL_KERNEL_CONSTANT_H

namespace CGAL{
  typedef float HQ_NT;//type to represent the index of one hquadrant
  enum Fct_type{TAN, COT, FIXED, TAG_M2};
  enum Circle_type{NORMAL,THREADED,POLAR,BIPOLAR};
  
  inline Fct_type auto_ftype(const HQ_NT& hquad){
  if (hquad >7 || hquad<2 || (hquad>3 && hquad<6))
    return TAN;
  return COT;
  };    
  
  #warning introduce this in a better manner
  template<class FT,class Point_3>
  inline FT compute_a(const Point_3& c,const FT& R2,const FT& squared_radius){
    return c.x() * c.x() + c.y() * c.y()+ c.z() * c.z() + R2 - squared_radius;
  };
  
  template<class FT>
  FT circle_center_coefficent(const FT& x,const FT& y,const FT& z,const FT& r2,const FT& R2){
    return ((FT)(0.5) +  (R2 - r2)/(FT)(2* (x*x +y*y +z*z))) ;
  }
  
  template<class SK, class circle_on_sphere>
  CGAL::Circle_type classify_one_circle(const circle_on_sphere& C){
    if (C.supporting_sphere().center().z()==0){
      typename SK::Point_3 Pt=C.center();
      if (Pt.z()==0 && Pt.y()==0 && Pt.x()==0)
      return CGAL::BIPOLAR;
    }
    std::vector<CGAL::Object> cont;
    typename SK::Plane_3 Pl=SK().construct_plane_3_object()(typename SK::Algebraic_kernel::Polynomial_1_3(0,1,0,0));
    typename SK::Intersect_3()(C.reference_sphere(),C.supporting_sphere(),Pl,std::back_inserter(cont));
    
    switch (cont.size()){
      case 0: 
        return CGAL::NORMAL;
      case 2:{ 
        std::pair<typename SK::Circular_arc_point_3,unsigned> p1,p2;
        CGAL::assign(p1,cont[0]);CGAL::assign(p2,cont[1]);
        CGAL::Sign s1=CGAL::sign(p1.first.x());
        CGAL::Sign s2=CGAL::sign(p2.first.x());
        if (s1==CGAL::opposite(s2))
          return CGAL::THREADED;
        else
          if (s1!=s2) return CGAL::POLAR;
        }
        break;
      
      case 1:
        if (CGAL::abs(C.extremal_point_z())==C.reference_sphere().radius()) return CGAL::POLAR;
    }
    return CGAL::NORMAL;
  }  
  
  
}

#endif
