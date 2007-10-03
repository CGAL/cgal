#ifndef CGAL_SPHERE_WITH_RADIUS_3_H
#define CGAL_SPHERE_WITH_RADIUS_3_H

namespace CGAL {
  template <class SK> 
    class Sphere_with_radius_3
    : public SK::Kernel_base::Sphere_with_radius_3
  {
    typedef typename SK::FT                    FT;
    typedef typename SK::Point_3               Point_3;
    typedef typename SK::Plane_3               Plane_3;
    typedef typename SK::Kernel_base::Sphere_with_radius_3 RSphere_with_radius_3; 
   
  
  public:
    typedef  RSphere_with_radius_3 Repd;
    typedef  SK   R; 
    
    const Repd& rep() const{
      return *this;
    }
    
    Repd& rep(){
      return *this;
    }
    
    Sphere_with_radius_3(const RSphere_with_radius_3& c)
    :RSphere_with_radius_3(c){}
    
    Sphere_with_radius_3()
      : RSphere_with_radius_3(typename R::Construct_sphere_with_radius_3()())
      {}

    Sphere_with_radius_3(const FT& _r,const Point_3& _c)
      : RSphere_with_radius_3(typename R::Construct_sphere_with_radius_3()(_r,_c))
      {}

    typename Qualified_result_of
    <typename R::Construct_radius_sphere_with_radius_3, Sphere_with_radius_3>::type
    //const Sphere_3 &
    radius() const
    {
      return typename R::Construct_radius_sphere_with_radius_3()(*this);
    }
  };
}
#endif //CGAL_SPHERE_WITH_RADIUS_3_H
