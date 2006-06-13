#ifndef CGAL_SPHERE_COORDINATE_H
#define CGAL_SPHERE_COORDINATE_H
#include <iostream>
struct Coordinate_index {
  Coordinate_index(): c_(-1){}
  explicit Coordinate_index(int c): c_(c){
    CGAL_assertion(c_ >=0);
  }
  /*int index() const {
    return c_;
    }*/
  operator unsigned int () const {
    CGAL_precondition(is_valid());
    return c_;
  }
  bool is_valid() const {
    return c_>=0 && c_ <3;
  }
  bool operator==(const Coordinate_index &o) const {
    return c_== o.c_;
  }
  bool operator!=(const Coordinate_index &o) const {
    return c_!= o.c_;
  }
  std::ostream &write(std::ostream &out) const {
    out << c_;
    return out;
  }
  int c_;
};

inline const Coordinate_index operator-(int one, Coordinate_index c) {
  CGAL_precondition(one==1);
  return Coordinate_index(one-static_cast<unsigned int>(c));
}

inline std::ostream &operator<<(std::ostream &out, Coordinate_index c) {
  return c.write(out);
}

#endif
