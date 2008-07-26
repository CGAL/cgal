#ifndef CGAL_PLANE_3_BBOX_3_DO_INTERSECT_H
#define CGAL_PLANE_3_BBOX_3_DO_INTERSECT_H

// Opcode like

CGAL_BEGIN_NAMESPACE

namespace CGALi {

  template <class K>
  bool do_intersect(const typename K::Plane_3& plane, 
    const CGAL::Bbox_3& bbox,
    const K& kernel)
  {	
    K::Point_3 p_max, p_min;
    get_min_max<K>(plane.orthogonal_vector(), bbox, p_min, p_max);
    return ! (plane.oriented_side(p_max) == ON_NEGATIVE_SIDE || 
      plane.oriented_side(p_min) == ON_POSITIVE_SIDE);
  }

  template <class K>
  void get_min_max(const typename K::Vector_3& p,
    const CGAL::Bbox_3& bbox,
    typename K::Point_3& p_min, 
    typename K::Point_3& p_max)
  {
    if(p.x() > 0) {
      if(p.y() > 0) {
	if(p.z() > 0) { p_min = K::Point_3(bbox.xmin(), bbox.ymin(),bbox.zmin()); 
	p_max = K::Point_3(bbox.xmax(), bbox.ymax(),bbox.zmax());}
	else {							     p_min = K::Point_3(bbox.xmin(), bbox.ymin(),bbox.zmax()); 
	p_max = K::Point_3(bbox.xmax(), bbox.ymax(),bbox.zmin());}
      }
      else {
	if(p.z() > 0) { p_min = K::Point_3(bbox.xmin(), bbox.ymax(),bbox.zmin()); 
	p_max = K::Point_3(bbox.xmax(), bbox.ymin(),bbox.zmax());}
	else {					         p_min = K::Point_3(bbox.xmin(), bbox.ymax(),bbox.zmax()); 
	p_max = K::Point_3(bbox.xmax(), bbox.ymin(),bbox.zmin());}
      }
    }
    else {
      if(p.y() > 0) {
	if(p.z() > 0) { p_min = K::Point_3(bbox.xmax(), bbox.ymin(),bbox.zmin()); 
	p_max = K::Point_3(bbox.xmin(), bbox.ymax(),bbox.zmax());}
	else {					         p_min = K::Point_3(bbox.xmax(), bbox.ymin(),bbox.zmax()); 
	p_max = K::Point_3(bbox.xmin(), bbox.ymax(),bbox.zmin());}
      }
      else {
	if(p.z() > 0) { p_min = K::Point_3(bbox.xmax(), bbox.ymax(),bbox.zmin()); 
	p_max = K::Point_3(bbox.xmin(), bbox.ymin(),bbox.zmax());}
	else {					         p_min = K::Point_3(bbox.xmax(), bbox.ymax(),bbox.zmax()); 
	p_max = K::Point_3(bbox.xmin(), bbox.ymin(),bbox.zmin());}
      }
    }
  }

} // namespace CGALi

template <class K>
bool do_intersect(const CGAL::Plane_3<K>& plane, 
		  const CGAL::Bbox_3& bbox)
{
  return typename K::Do_intersect_3()(plane, bbox);
}

template <class K>
bool do_intersect(const CGAL::Bbox_3& bbox, 
		  const CGAL::Plane_3<K>& plane)
{
  return typename K::Do_intersect_3()(plane, bbox);
}


CGAL_END_NAMESPACE

#endif  // CGAL_PLANE_3_BBOX_3_DO_INTERSECT_H


