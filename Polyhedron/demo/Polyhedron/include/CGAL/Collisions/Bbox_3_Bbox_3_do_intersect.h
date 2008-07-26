#ifndef CGAL_BBOX_3_BBOX_3_DO_INTERSECT_H
#define CGAL_BBOX_3_BBOX_3_DO_INTERSECT_H

//#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

#undef min
#undef max

namespace CGALi {

  // assumes that the intersection with the supporting plane has
  // already been checked.
  template <class K>
  bool do_intersect(const CGAL::Bbox_3& c1, 
    const CGAL::Bbox_3& c2,
    const K& kernel)
  {
    for(int i = 0; i < 3; ++i) {
      if(c1.max(i) < c2.min(i) || c1.min(i) > c2.max(i))
	return false;
    }
    return true;
  }

} // namespace CGALi

template <class K>
bool do_intersect(const CGAL::Bbox_3& c, 
		  const CGAL::Bbox_3& bbox)
{
  return typename K::Do_intersect_3()(c, bbox);
}

CGAL_END_NAMESPACE

#endif  // CGAL_BBOX_3_BBOX_3_DO_INTERSECT_H
