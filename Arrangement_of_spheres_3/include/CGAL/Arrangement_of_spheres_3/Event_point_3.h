#ifndef CGAL_ARRANGEMENT_OF_SPHERES_EP_H
#define CGAL_ARRANGEMENT_OF_SPHERES_EP_H
#include <CGAL/Arrangement_of_spheres_3/Sphere_line_intersection.h>
#include <CGAL/Kinetic/basic.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE 

template <class Tr>
class Event_point_3: public Sphere_line_intersection<Tr> {
  typedef Event_point_3 This;
  typedef Sphere_line_intersection<Tr> P;
public:
  Event_point_3(){}
  Event_point_3(P sp3): P(sp3) {CGAL_assertion(sp3.is_valid());}
  Event_point_3(typename P::NT t): P(typename P::Sphere_3(sweep_point<typename P::Point_3>(t), 0), 
				     typename P::Line_3(typename P::Point_3(0,0,0), 
							sweep_vector<typename P::Vector_3>())){

    //std::cout << P::sphere() << std::endl;
    //std::cout << P::line() << std::endl;
  }

  bool operator==(const This &o) const {
    return P::compare(o, sweep_coordinate()) == CGAL::EQUAL;
  }
  bool operator!=(const This &o) const {
    return P::compare(o, sweep_coordinate()) != CGAL::EQUAL;
  }
  bool operator<=(const This &o) const {
    return P::compare(o, sweep_coordinate()) != CGAL::LARGER;
  }
  bool operator>=(const This &o) const {
    return P::compare(o, sweep_coordinate()) != CGAL::SMALLER;
  }
  bool operator<(const This &o) const {
    return P::compare(o, sweep_coordinate()) == CGAL::SMALLER;
  }
  bool operator>(const This &o) const {
    return P::compare(o, sweep_coordinate()) == CGAL::LARGER;
  }

  This operator-() const {
    return This(this->flip_on_sweep());
  }


  std::pair<double, double> interval_value() const {
    return P::interval_coordinate(sweep_coordinate());
  }

  double double_value() const {
    return P::approximate_coordinate(sweep_coordinate());
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

#endif
