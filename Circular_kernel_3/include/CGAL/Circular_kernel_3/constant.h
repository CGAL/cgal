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
    return ((FT)(0.5) +  (R2 - r2)/(FT)(2* x*x +y*y +z*z)) ;    
  }
  
}

#endif
