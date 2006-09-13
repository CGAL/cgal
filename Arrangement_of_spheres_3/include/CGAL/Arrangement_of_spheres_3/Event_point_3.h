#ifndef CGAL_ARRANGEMENT_OF_SPHERES_H
#define CGAL_ARRANGEMENT_OF_SPHERES_H
#include <CGAL/Arrangement_of_spheres_3/Sphere_line_intersection.h>
#include <CGAL/Arrangement_of_spheres_3/coordinates.h>

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

namespace CGAL {
  template <class Tr>
  std::pair<double,double> to_interval(const ::Event_point_3<Tr> &a) {
    return a.interval_value();
  }
  template <class Tr>
  double to_double(const ::Event_point_3<Tr> &a) {
    return a.double_value();
  }
};

namespace std {
  
  template <class Tr>
  struct numeric_limits<Event_point_3<Tr> > {
    typedef Event_point_3<Tr> T;
    static const bool is_specialized = true;
    static T min() throw () {CGAL_assertion(0); return T(0);}
    static T max() throw () {CGAL_assertion(0); return T(0);}
    static const int digits =0;
    static const int digits10 =0;
    static const bool is_signed = true;
    static const bool is_integer = false;
    static const bool is_exact = true;
    static const int radix =0;
    static T epsilon() throw(){return T(0);}
    static T round_error() throw(){return T(0);}
    static const int min_exponent=0;
    static const int min_exponent10=0;
    static const int max_exponent=0;
    static const int max_exponent10=0;
    static const bool has_infinity=false;
    static const bool has_quiet_NaN = true;
    static const bool has_signaling_NaN= false;
    static const float_denorm_style has_denorm= denorm_absent;
    static const bool has_denorm_loss = false;
    static T infinity() throw() {CGAL_assertion(0); return T(0);}
    static T quiet_NaN() throw(){return T();}
    static T denorm_min() throw() {return T(0);}
    static const bool is_iec559=false;
    static const bool is_bounded =false;
    static const bool is_modulo= false;
    static const bool traps = false;
    static const bool tinyness_before =false;
    static const float_round_style round_stype = round_toward_zero;
  };

}

#endif
