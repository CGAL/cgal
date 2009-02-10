#ifndef CGAL_AOS3_CENTER_POINT_3_H
#define CGAL_AOS3_CENTER_POINT_3_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>


CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

template <class K, class C>
struct Center_point_3 {
  typedef Center_point_3<K,C> This;
  typedef K Key;
  typedef C Coord;
  Center_point_3(K k,
		 const Coord &pt): k_(k), t_(pt){}
  Key key() const {return k_;}
  const Coord &coord() const {return t_;}
  
  std::ostream &write(std::ostream&out) const {
    return out << k_ << ": " << t_ << std::endl;
  }

  CGAL_COMPARISONS2(k_, t_);

private:
  Key k_;
  Coord t_;
};


CGAL_OUTPUT2(Center_point_3);
CGAL_AOS3_END_INTERNAL_NAMESPACE
#endif
