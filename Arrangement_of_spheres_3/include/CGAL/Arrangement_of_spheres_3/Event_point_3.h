#ifndef CGAL_ARRANGEMENT_OF_SPHERES_EP_H
#define CGAL_ARRANGEMENT_OF_SPHERES_EP_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Tools/utility_macros.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE 

template <class Tr>
class Event_point_3: public Tr::Sphere_point_3 {
  typedef Event_point_3 This;
  typedef typename Tr::Sphere_point_3 P;
public:
  Event_point_3(){}
  Event_point_3(P sp3): P(sp3) {CGAL_assertion(sp3.is_valid());}
  Event_point_3(typename P::NT t): P(sweep_point<typename P::Point_3>(t)){

    //std::cout << P::sphere() << std::endl;
    //std::cout << P::line() << std::endl;
  }
  
  CGAL::Comparison_result compare(const This &o) const {
    return P::compare_c(o, sweep_coordinate());
  }
  
  //using P::compare;

  CGAL_COMPARISONS;
  
 
  This operator-() const {
    return This(this->flip_on_sweep());
  }


  std::pair<double, double> approximating_interval(double=0) const {
    CGAL_precondition(sweep_coordinate().index()==0);
    return std::make_pair(P::bbox().xmin(), P::bbox().xmax());
  }

  double approximation(double =0) const {
    return .5*(approximating_interval().first+approximating_interval().second);
  }
  
};

/*template <class Tr>
std::pair<double,double> to_interval(const ::Event_point_3<Tr> &a) {
  return a.interval_value();
}

template <class Tr>
double to_double(const ::Event_point_3<Tr> &a) {
  return a.double_value();
  }*/


CGAL_AOS3_END_INTERNAL_NAMESPACE
CGAL_REAL_EMBEDDABLE1(CGAL_AOS3_INTERNAL_NS::Event_point_3);
//CGAL_COMPARABLE1(CGAL_AOS3_INTERNAL_NS::Event_point_3);

#endif
