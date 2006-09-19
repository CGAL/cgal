#include <CGAL/Arrangement_of_spheres_3/Slice.h>

#define DPRINT(x) x


/* 
   predictes--------------------------------------------------------
*/






int Slice::sphere_location(const T::Sphere_point_3& sp,
			   T::Key s) const {
  int r=0;
  CGAL::Bounded_side bs=t_.bounded_side_of_sphere(s,
						  sp);
							      
  /*T::Event_point_3 ep=sphere_start(locate_point);
  T::Line_3 l= ep.line();
  T::Event_point_3 a(sphere(s), l);
  if (a.is_valid()) {
    T::Event_point_3 b(sphere(s), l.opposite());
    CGAL::Comparison_result ca= a.compare(ep, 2);
    CGAL::Comparison_result cb= b.compare(ep, 2);*/
	
  if (bs== CGAL::ON_BOUNDED_SIDE) {
    r |= Sds::Curve::lIN_BIT;
  } else if (bs==CGAL::ON_BOUNDARY) {
    r |= Sds::Curve::lIN_BIT;
    r |= Sds::Curve::lOUT_BIT;
  } else {
    r |= Sds::Curve::lOUT_BIT;
  }

  /*
    NOTE can make it more specific for intersection point location
  */

  CGAL::Comparison_result xo= t_.compare_sphere_center_c(s, sp,
							  plane_coordinate(0));
  CGAL::Comparison_result yo= t_.compare_sphere_center_c(s, sp, 
							 plane_coordinate(1));
  if (xo  != CGAL::LARGER) {
    r |= Sds::Curve::lR_BIT;
  } 
  if (xo != CGAL::SMALLER){
    r |= Sds::Curve::lL_BIT;
  }
  if (yo != CGAL::LARGER) {
    r |= Sds::Curve::lT_BIT;
  } 
  if (yo != CGAL::SMALLER) {
    r |= Sds::Curve::lB_BIT;
  }
  return r;
}




 




bool Slice::behind_arc(const T::Sphere_point_3 &ep,
		       Sds::Curve arc,
		       int location) const{
  CGAL_precondition(!arc.is_inside());
  CGAL_precondition(arc.is_compatible_location(location));
  T::Coordinate_index C=arc.is_weakly_incompatible(location);

  if (C.is_valid()) { // we know it is opposite the center
    /*CGAL::Comparison_result cr = t_.compare_sphere_center_c(arc.key(),
							    ep, C);
    if (cr == CGAL::SMALLER && arc.is_negative()
    || cr == CGAL::LARGER && !arc.is_negative()) return false;*/
    bool b= t_.is_over_circle_c(arc.key(), ep, C);
    return b;
  } else {
    return false;
  }
}


void Slice::point_sphere_orientation(const T::Sphere_point_3 &time,
				     T::Key sphere,
				     std::vector<int> &locations
				     /*,
				       std::vector<Sds::Curve> &edges*/) const {
  if (locations[sphere.input_index()]==0){
    locations[sphere.input_index()]=sphere_location(time, sphere);
    //DPRINT(std::cout << "For sphere " << sphere << " got " << SL::decode(locations[sphere]) << std::endl);
  }
}
