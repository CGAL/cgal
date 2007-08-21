#ifndef CGAL_SPHERE_COORDINATE_H
#define CGAL_SPHERE_COORDINATE_H
//#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <iostream>


CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

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
  /*Coordinate_index(Coordinate_X):c_(0){}
  Coordinate_index(Coordinate_Y):c_(1){}
  Coordinate_index(Coordinate_Z):c_(2){}*/
  
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


class Coordinate_X;
class Coordinate_Y;
class Coordinate_Z;

struct Coordinate_X {
  static int index() {return 0;}
  static Coordinate_index object() {return Coordinate_index(0);}
  template <class V>
  static const V& basis_3() {
    static V v(1,0,0);
    return v;
  }
};

struct Coordinate_Y {
  static int index() {return 1;}
  typedef Coordinate_Z Other_plane_coordinate;
  static int plane_index() {return 0;}
  static Coordinate_index object() {return Coordinate_index(1);}
  template <class V>
  static const V& basis_3() {
    static V v(0,1,0);
    return v;
  }
};

struct Coordinate_Z {
  static int index() {return 2;}
  typedef Coordinate_Y Other_plane_coordinate;
  static int plane_index() {return 1;}
  static Coordinate_index object() {return Coordinate_index(2);}
  template <class V>
  static const V& basis_3() {
    static V v(0,0,1);
    return v;
  }
};


typedef Coordinate_X Sweep_coordinate;
typedef Coordinate_Y Plane_coordinate_0;
typedef Coordinate_Z Plane_coordinate_1;


inline Coordinate_index sweep_coordinate() {
  return Sweep_coordinate::object();
}
inline Coordinate_index plane_coordinate(int i) {
  if (i ==0) return Plane_coordinate_0::object();
  else return Plane_coordinate_1::object();
}
inline unsigned int project(Coordinate_index i) {
  if (i == Plane_coordinate_0::object()) return 0;
  else {
    CGAL_assertion(i== Plane_coordinate_1::object());
    return 1;
  }
}

inline Coordinate_index other_plane_coordinate(Coordinate_index i) {
  if (i == Plane_coordinate_0::object()) return Plane_coordinate_1::object();
  else {
    CGAL_assertion(i== Plane_coordinate_1::object());
    return Plane_coordinate_0::object();
  }
}

template <class Vector>
inline Vector sweep_vector() {
  int v[3]={0,0,0};
  v[sweep_coordinate().index()]=1;
  return Vector(v[0], v[1], v[2]);
}

template <class Point, class NT>
inline Point sweep_point(NT n) {
  NT v[3]={0,0,0};
  v[sweep_coordinate().index()]=n;
  return Point(v[0], v[1], v[2]);
}

template <class K>
class Unprojector {
  typename K::FT z_;
public:
  Unprojector(typename K::FT z): z_(z){}
 
  typename K::Point_3 operator()(const typename K::Point_2 &pt) const {
    typename K::FT v[3];
    v[Sweep_coordinate::index()]= z_;
    v[Plane_coordinate_0::index()]= pt.x();
    v[Plane_coordinate_1::index()]= pt.y();
    return typename K::Point_3(v[0], v[1], v[2]);
  }

  typename K::Vector_3 operator()(const typename K::Vector_2 &pt) const {
      typename K::FT v[3];
      v[Sweep_coordinate::index()]= 0;
      v[Plane_coordinate_0::index()]= pt.x();
      v[Plane_coordinate_1::index()]= pt.y();
      return typename K::Vector_3(v[0], v[1], v[2]);
  }
  typename K::Point_2 operator()(const typename K::Point_3 &pt) const {
    return typename K::Point_2(pt[Plane_coordinate_0::index()],
			       pt[Plane_coordinate_1::index()]);
  }

  typename K::Vector_2 operator()(const typename K::Vector_3 &pt) const {
     return typename K::Vector_2(pt[Plane_coordinate_0::index()],
				 pt[Plane_coordinate_1::index()]);
  }
  
  Bbox_2 operator()(const Bbox_3 &bb) const {
    typename K::Point_3 minp(bb.xmin(), bb.ymin(), bb.zmin());
    typename K::Point_3 maxp(bb.xmax(), bb.ymax(), bb.zmax());
    typename K::Point_2 minp2=operator()(minp);
    typename K::Point_2 maxp2=operator()(maxp);
    return minp2.bbox()+ maxp2.bbox();
  }
};



CGAL_AOS3_END_INTERNAL_NAMESPACE
#endif
