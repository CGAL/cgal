#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Coordinate_index.h>
CGAL_BEGIN_NAMESPACE

CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_sphere_center_c(Sphere_3_key a,
									const Sphere_point_3& d,
									Coordinate_index C) const {
  return CGAL::Comparison_result(-d.compare(table_->center(a), C));
}

CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_sphere_center_c(Sphere_3_key a,
									FT d,
									Coordinate_index C) const {
  return CGAL::compare(table_->center(a)[C.index()], d);
}

CGAL_AOS3_TEMPLATE
CGAL::Comparison_result
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_sphere_centers_c(Sphere_3_key a, Sphere_3_key b,
									 Coordinate_index C) const {
  return CGAL::compare(table_->center(a)[C.index()], table_->center(b)[C.index()]);
}



CGAL_AOS3_TEMPLATE
bool 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::sphere_intersects_sweep(Sphere_3_key sphere,
									const Sphere_point_3& ep) const {
  Event_pair be= sphere_events(sphere);
  if (be.first.compare(ep, CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) != CGAL::LARGER 
      && be.second.compare(ep, CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) != CGAL::SMALLER) return true;
  else return false;
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_pair 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::sphere_events(Sphere_3_key ind) const {
  FT v[3];
  v[CGAL_AOS3_INTERNAL_NS::plane_coordinate(0).index()]=0;
  v[CGAL_AOS3_INTERNAL_NS::plane_coordinate(1).index()]=0;
  v[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()]=1;
  return Event_pair(Sphere_point_3(table_->sphere(ind), 
				   Line_3(table_->center(ind), Vector_3(v[0], v[1], v[2]))),
		    Sphere_point_3(table_->sphere(ind), 
				   Line_3(table_->center(ind), Vector_3(-v[0], -v[1], -v[2]))));
}


CGAL_END_NAMESPACE
