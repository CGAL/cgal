#include <CGAL/Arrangement_of_spheres_3/Slice.h>

#define DPRINT(x) x


/* 
   predictes--------------------------------------------------------
*/






int Slice::sphere_location(const T::Sphere_point_3& sp,
			   T::Key locate_point,
			   T::Key s) const {
  int r=0;
  CGAL::Bounded_side bs=t_.bounded_side_of_sphere(s,
						  locate_point,
						  locate_point,
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
  
  /*T::Plane_3 lrp(center(s), T::Vector_3(1, 0, 0));
  T::Plane_3 tbp(center(s), T::Vector_3(0, 1, 0));
  CGAL::Comparison_result xo= oriented_side(lrp, ep);
  CGAL::Oriented_side yo= oriented_side(tbp, ep);*/
  //CGAL::Comparison_result xo= compare_event_point_to_rule(ep, s, 0);
  //CGAL::Comparison_result yo= compare_event_point_to_rule(ep, s, 1);
  CGAL::Comparison_result xo= t_.compare_sphere_centers_c(s, locate_point,
							  plane_coordinate(0));
  CGAL::Comparison_result yo= t_.compare_sphere_centers_c(s, locate_point, 
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
		       T::Key ind, Sds::Curve arc,
		       int location) const{
  CGAL_assertion(arc.is_compatible_location(location));
  T::Coordinate_index C=arc.is_weakly_incompatible(location);

  if (C.is_valid()) {
    CGAL::Bounded_side bs;
    if (C== plane_coordinate(0)) {
      bs= t_.bounded_side_of_sphere(arc.key(),
				    arc.key(),
				    ind,
				    ep);
    } else {
      bs= t_.bounded_side_of_sphere(arc.key(),
				    ind,
				    arc.key(),
				    ep);
    }
    /*NT v[2];
    v[C]= t_.center_c(arc.key(), C);
    v[1-C]= t_.center_c(ind, 1-C);
    T::Line_3 l(T::Point_3(v[0], v[1], 0), 
		T::Vector_3(0,0,1));
    T::Event_point_3 fp(sphere(arc.key()),
			l);
    if (!fp.is_valid()) return false;
    T::Event_point_3 bp(sphere(arc.key()),
			l.opposite());
    if (fp <= ep && bp >=ep) {
      return true;
    } else {return false;}
    } else return false;*/
    return bs== CGAL::ON_BOUNDED_SIDE;
  } else {
    return false;
  }
}


void Slice::point_sphere_orientation(const T::Sphere_point_3 &time,
				     T::Key front_point,
				     T::Key sphere,
				     std::vector<int> &locations
				     /*,
				       std::vector<Sds::Curve> &edges*/) const {
  if (locations[sphere.input_index()]==0){
    locations[sphere.input_index()]=sphere_location(time, front_point, sphere);
    //DPRINT(std::cout << "For sphere " << sphere << " got " << SL::decode(locations[sphere]) << std::endl);
  }
}
