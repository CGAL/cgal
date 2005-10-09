#ifndef CGAL_KDS_INTERNAL_DELAUNAY_LIFTING_H
#define CGAL_KDS_INTERNAL_DELAUNAY_LIFTING_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE

template <class K>
struct Delaunay_lifting {
  typedef typename K::Point_3 argument_type;
  typedef typename K::Motion_function result_type;
  
  result_type operator()(const argument_type &p) const {
    return CGAL::square(p.x())+ CGAL::square(p.y()) + CGAL::square(p.z());
  }
  result_type operator()(const typename K::Weighted_point_3 &wp) const {
    return operator()(wp.point()) - wp.weight();
  }
};

CGAL_KDS_END_INTERNAL_NAMESPACE

#endif
