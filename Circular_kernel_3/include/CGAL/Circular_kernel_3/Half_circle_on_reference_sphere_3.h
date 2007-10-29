#ifndef CGAL_SPHERICAL_HALF_CIRCLE_ON_REFERENCE_SPHERE_H
#define CGAL_SPHERICAL_HALF_CIRCLE_ON_REFERENCE_SPHERE_H

namespace CGAL {
  namespace CGALi {

  template <class SK>
  class Half_circle_on_reference_sphere_3
  {
    protected:
    const typename SK::Circle_on_reference_sphere_3& C;
    CGAL::Hcircle_type pos;
    public:
    const typename SK::Circle_on_reference_sphere_3& supporting_circle() const {return C;}
    const CGAL::Hcircle_type& get_position() const {return pos;}
    Half_circle_on_reference_sphere_3(const typename SK::Circle_on_reference_sphere_3& C,CGAL::Hcircle_type pos):C(C),pos(pos){}
  };
    
  }
}

#endif
