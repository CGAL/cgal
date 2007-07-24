#ifndef CGAL_SPHERE_COORDINATE_H
#define CGAL_SPHERE_COORDINATE_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <iostream>


CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

class Coordinate_X;
class Coordinate_Y;
class Coordinate_Z;

struct Coordinate_X {
  static int index() {return 0;}
  typedef Coordinate_Y Other_plane_coordinate;
};

struct Coordinate_Y {
  static int index() {return 1;}
  typedef Coordinate_X Other_plane_coordinate;
};

struct Coordinate_Z {
  static int index() {return 2;}
};


class Coordinate_index {
  int c_;
public:
  typedef Coordinate_index This;
  Coordinate_index(): c_(-1){}
  //template <class C>
  explicit Coordinate_index(int c): c_(c){
    CGAL_assertion(c_ >=0 && c_ < 3);
    // bool check_constructor;
  }
  Coordinate_index(Coordinate_X):c_(0){}
  Coordinate_index(Coordinate_Y):c_(1){}
  Coordinate_index(Coordinate_Z):c_(2){}
  
  /*explicit Coordinate_index(int c, bool): c_(c){
    CGAL_assertion(c_ >=0 && c_<3);
    }*/
  /*int index() const {
    return c_;
    }*/
  /*operator unsigned int () const {
    CGAL_precondition(is_valid());
    return c_;
    }*/
  CGAL_GETNR(unsigned int, index,  CGAL_precondition(is_valid());
			return c_);
  CGAL_IS(valid,  return c_>=0 && c_ <3);

  CGAL_COMPARISONS1(c_);

  /*bool operator==(const Coordinate_index &o) const {
    return c_== o.c_;
  }
  bool operator!=(const Coordinate_index &o) const {
    return c_!= o.c_;
    }*/

  std::ostream &write(std::ostream &out) const {
    char c[]="XYZ";
    out << c[c_];
    return out;
  }
  static Coordinate_index X() {
     Coordinate_index r;
     r.c_= 0;
     return r;
  }
  static Coordinate_index Y() {
    Coordinate_index r;
     r.c_= 1;
     return r;
  }
  static Coordinate_index Z() {
    Coordinate_index r;
    r.c_= 2;
    return r;
  }

};

/*inline const Coordinate_index operator-(int one, Coordinate_index c) {
  CGAL_precondition(one==1);
  return Coordinate_index(one-static_cast<unsigned int>(c));
  }*/

CGAL_OUTPUT(Coordinate_index);



CGAL_AOS3_END_INTERNAL_NAMESPACE
#endif
