#include <CGAL/Arrangement_of_spheres_3/Slice.h>

#define DPRINT(x) x


/* 
   predictes--------------------------------------------------------
*/

bool Slice::intersects_rule(int s, int rule, int C) const {
  NT d= spheres_[s].center()[C]-spheres_[rule].center()[C];
  if (CGAL::square(d) <= spheres_[s].squared_radius()) return true;
  else return false;
}

CGAL::Comparison_result Slice::compare_equipower_point_to_rule(int a, int b,
							       int c, int C) const{
  // NOTE redo this with computing it directly
  T::Plane_3 eqp= equipower_plane(a,b);
  T::Line_3 l(center(a), (center(a)-center(b)));
  CGAL::Object o= CGAL::intersection(eqp, l);
  T::Point_3 pt;
  CGAL_assertion(CGAL::assign(pt, o));
  CGAL::assign(pt, o);
  return CGAL::compare(pt[C], center_c(c, C));

}


CGAL::Comparison_result  Slice::compare_sphere_center_to_rule(int sphere, 
							      int rule_sphere,
							      int C) const {
  return CGAL::compare(spheres_[sphere].center()[C],
		       spheres_[rule_sphere].center()[C]);
}


CGAL::Sign Slice::sign_of_separating_plane_normal_c(int a, int b, int C) const {
  CGAL::Sign sn;
  if (C==0) {
    sn = CGAL::sign(-spheres_[b].center()[1] + spheres_[a].center()[1]);
  } else {
    sn = CGAL::sign( spheres_[b].center()[0] - spheres_[a].center()[0]);
  }
  return sn;
}




CGAL::Sign Slice::sign_of_equipower_plane_normal_c(int a, 
						   int b, int C) const {
  CGAL::Sign sn= CGAL::sign(spheres_[a].center()[C]-spheres_[b].center()[C]);
  CGAL_assertion(sn== CGAL::sign(equipower_plane(a,b).orthogonal_vector()[C]));
  return sn;
}



CGAL::Oriented_side Slice::oriented_side_of_equipower_plane(int a, int b,
							    const T::Sphere_point_3 &s) const {
  CGAL_assertion(s.is_valid());
  T::Plane_3 p= equipower_plane(a,b);
  CGAL::Object o= CGAL::intersection(p, s.line());
  T::Point_3 pt;
  T::Line_3 l;
  if (CGAL::assign(pt, o)){
    CGAL::Comparison_result c= s.compare_on_line(pt);
    NT d= p.orthogonal_vector() * s.line().to_vector();
    if (d > 0) return CGAL::enum_cast<CGAL::Oriented_side>(c);
    else return CGAL::enum_cast<CGAL::Oriented_side>(-c);
  } else if (CGAL::assign( l, o)){
    return CGAL::ON_ORIENTED_BOUNDARY;
  } else {
    return p.oriented_side(s.line().point());
  }
}


CGAL::Oriented_side Slice::oriented_side_of_center_plane(int a, int b,
							 int sphere_center) const {
  T::Vector_3 d(spheres_[b].center()-spheres_[a].center());
  T::Line_2 l(T::Point_2(spheres_[a].center().x(),
			 spheres_[a].center().y()), 
	      T::Vector_2(d.x(), d.y()));
  return l.oriented_side(T::Point_2(spheres_[sphere_center].center().x(),
				    spheres_[sphere_center].center().y()));
}

CGAL::Comparison_result Slice::compare_sphere_centers_c(int a, int b, int C) const {
  return CGAL::compare(spheres_[a].center()[C], spheres_[b].center()[C]);
}





















int Slice::sphere_location(int locate_point, int s) const {
  int r=0;
  T::Event_point_3 ep=sphere_start(locate_point);
  T::Line_3 l= ep.line();
  T::Event_point_3 a(sphere(s), l);
  if (a.is_valid()) {
    T::Event_point_3 b(sphere(s), l.opposite());
    CGAL::Comparison_result ca= a.compare(ep, 2);
    CGAL::Comparison_result cb= b.compare(ep, 2);
	
    if (ca == CGAL::SMALLER && cb == CGAL::LARGER) {
      r |= Sds::Curve::lIN_BIT;
    } else if (ca == CGAL::EQUAL || cb == CGAL::EQUAL) {
      r |= Sds::Curve::lIN_BIT;
      r |= Sds::Curve::lOUT_BIT;
    } else {
      r |= Sds::Curve::lOUT_BIT;
    }
    //std::cout << ca << " " << cb << std::endl;
  } else {
    r |= Sds::Curve::lOUT_BIT;
  }
  
  T::Plane_3 lrp(center(s), T::Vector_3(1, 0, 0));
  T::Plane_3 tbp(center(s), T::Vector_3(0, 1, 0));
  CGAL::Oriented_side xo= oriented_side(lrp, ep);
  CGAL::Oriented_side yo= oriented_side(tbp, ep);
  if (compare_sphere_center_to_rule BLAHHHHHHHHH != CGAL::SMALLER) {
    r |= R_BIT;
  } 
  if (xo != CGAL::LARGER){
    r |= L_BIT;
  }
  if (yo != CGAL::SMALLER) {
    r |= T_BIT;
  } 
  if (yo != CGAL::LARGER) {
    r |= B_BIT;
  }
  return r;
}




 




bool Slice::behind_arc(T::Event_point_3 ep, int ind, Sds::Curve arc,
		       int location) const{
  CGAL_assertion(arc.is_compatible_location(location));
  int C=arc.is_weakly_incompatible(location);
  if (C != -1) {
    NT v[2];
    v[C]= spheres_[arc.index()].center()[C];
    v[1-C]= spheres_[ind].center()[1-C];
    T::Line_3 l(T::Point_3(v[0], v[1], 0), 
		T::Vector_3(0,0,1));
    T::Event_point_3 fp(spheres_[arc.index()],
			l);
    if (!fp.is_valid()) return false;
    T::Event_point_3 bp(spheres_[arc.index()],
			l.opposite());
    if (fp <= ep && bp >=ep) {
      return true;
    } else {return false;}
  } else return false;
}


void Slice::point_sphere_orientation(int front_point,
				     int sphere,
				     std::vector<int> &locations
				     /*,
				       std::vector<Sds::Curve> &edges*/) const {
  typedef T::Sphere_location SL;
  if (locations[sphere]==0){
    locations[sphere]=sphere_location(front_point, sphere);
    //DPRINT(std::cout << "For sphere " << sphere << " got " << SL::decode(locations[sphere]) << std::endl);
  }
}
