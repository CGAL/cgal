#ifndef CGAL_SPHERICAL_KERNEL_SPHERE_WITH_RADIUS_3_H
#define CGAL_SPHERICAL_KERNEL_SPHERE_WITH_RADIUS_3_H

//~ #include <CGAL/Cartesian/Sphere_3.h>
#include <CGAL/Sphere_3.h>

namespace CGAL {
  namespace CGALi {
    
  template<class SK>
  class Sphere_with_radius_3:public CGAL::Sphere_3<SK>{
    typedef typename SK::Algebraic_kernel::Root_of_2 Root_of_2;
    typedef typename SK::FT FT;
    typedef typename SK::Algebraic_kernel AK;
    typedef typename SK::Sphere_3 Sphere_3;
    //---------------
    typename SK::template Handle<FT>::type  hrad;
    public:
    typedef typename SK::Point_3 Point_3;
    Sphere_with_radius_3():Sphere_3(){};
    Sphere_with_radius_3(const FT& _r,const Point_3& _c):Sphere_3(_c,_r*_r),hrad(_r){};
    const FT& radius() const {return CGAL::get(hrad);}
  };
    
  } // namespace CGALi
} // namespace CGAL

#endif //CGAL_SPHERICAL_KERNEL_SPHERE_WITH_RADIUS_3_H
