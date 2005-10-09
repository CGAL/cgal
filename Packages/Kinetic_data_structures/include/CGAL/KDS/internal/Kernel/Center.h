#ifndef CGAL_KDS_INTERNAL_CENTER_H
#define CGAL_KDS_INTERNAL_CENTER_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE

template <class K>
struct Center {
  typedef typename K::Point_3 argument_type;
  typedef typename K::Point_3 result_type;
  
  const result_type& operator()(const argument_type &p) const {
    return p;
  }
  const result_type& operator()(const typename K::Weighted_point_3 &wp) const {
    return wp.point();
  }
};

CGAL_KDS_END_INTERNAL_NAMESPACE

#endif
