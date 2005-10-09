#ifndef CGAL_KDS_CARTESIAN_MOVING_ORTHOSPHERE_H
#define CGAL_KDS_CARTESIAN_MOVING_ORTHOSPHERE_H
#include <CGAL/KDS/basic.h>
#include <CGAL/determinant.h>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE

template <class K>
class Cartesian_moving_orthosphere_3 {
  typedef typename K::Motion_function MF;
  typedef typename K::Point_3 PT;
public:
  template <class P>
  Cartesian_moving_orthosphere_3(const P &a,
				 const P &b,
				 const P &c,
				 const P &d, const K &k) {
    typename K::Delaunay_lifting lift= k.Delaunay_lifting_3_object();
    typename K::Center center k.center_3_object();
    initialize(center(a), lift(a),
	       center(b), lift(b),
	       center(c), lift(c),
	       center(d), lift(d));
  }
  template <class P>
  Cartesian_moving_orthosphere_3(const P &a,
				 const P &b,
				 const P &c, const K &k) {
    typename K::Delaunay_lifting lift= k.Delaunay_lifting_3_object();
    typename K::Center center k.center_3_object();
    initialize(center(a), lift(a),
	       center(b), lift(b),
	       center(c), lift(c));
  }
  
  const PT& center() const {
    return center_;
  }
  const MF l() const {
    return l_;
  };
  const MF o() const {
    return o_;
  };
  const PT &origin() 
protected:
  PT center_;
  PT origin_;
  MF l_, o_;
};

CGAL_KDS_END_INTERNAL_NAMESPACE

#endif
