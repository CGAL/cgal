#ifndef CGAL_SPHERICAL_CIRCLE_ON_REFERENCE_SPHERE_H
#define CGAL_SPHERICAL_CIRCLE_ON_REFERENCE_SPHERE_H

#include <CGAL/Circular_kernel_3/Circle_3.h>
#include <CGAL/Sphere_with_radius_3.h>
namespace CGAL {
  namespace CGALi {


    template<class SK>
    class Circle_on_reference_sphere_3
      : public CGAL::CGALi::Circle_3<SK,CGAL::CGALi::Circle_representation_3
          <typename CGAL::Sphere_with_radius_3<SK>,typename CGAL::Sphere_with_radius_3<SK>,SK > >{
    protected:
    typedef typename SK::Point_3 Point_3;
    typedef typename SK::FT FT;
    typedef typename SK::Sphere_with_radius_3 Sphere_3;            
    public:
      typedef CGAL::CGALi::Circle_3<SK,CGAL::CGALi::Circle_representation_3<
          typename CGAL::Sphere_with_radius_3<SK>,
          typename CGAL::Sphere_with_radius_3<SK>,SK > > Circle_3;
    protected:
      Circle_type _nature;//NORMAL,THREADED,POLAR,BIPOLAR

    public:
      
      Circle_on_reference_sphere_3(const FT& _r,const Point_3& _c,const Sphere_3& ref):Circle_3(Sphere_3(_r,_c),ref){
        _nature=SK::classify_one_circle(*this);
        #warning TRY WITH FLAG ONLY_SK_ALGEBRA
      }
      
      Circle_on_reference_sphere_3(const FT& _r,const Point_3& _c,Circle_type nat,const Sphere_3& ref):Circle_3(Sphere_3(_r,_c),ref),_nature(nat){}
        
      const Circle_type& type_of_circle_on_reference_sphere() const {  return _nature;}
      const FT& supporting_sphere_radius() const {  return this->supporting_sphere().radius();}
      const FT& supporting_sphere_squared_radius() const {  return this->supporting_sphere().squared_radius();}
      const Point_3& supporting_sphere_center() const {  return this->supporting_sphere().center();}    
      const Sphere_3& reference_sphere() const{return this->base.reference_sphere();}
      
      
      
      //according center of the circle=reference_sphere().center() + circle_center_coefficient() * (supporting_sphere().center()-reference_sphere().center())
      FT circle_center_coefficient() const {
        Point_3 center=this->supporting_sphere().center();
        return CGAL::circle_center_coefficent(center.x(),center.y(),center.z(),this->supporting_sphere().squared_radius(),reference_sphere().squared_radius());
      }
      
      FT extremal_point_z() const { return 2 / CGAL::compute_a<FT>(this->supporting_sphere().center(),
            reference_sphere().squared_radius(),this->supporting_sphere().squared_radius())
            * this->supporting_sphere().center().z() * reference_sphere().squared_radius();}

    };    
    
    
  }
}

#endif
