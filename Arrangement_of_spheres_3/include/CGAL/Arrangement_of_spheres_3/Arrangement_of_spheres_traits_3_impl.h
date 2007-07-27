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

CGAL_AOS3_TEMPLATE
 CGAL::Bounded_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::bounded_side_of_sphere(Sphere_3_key s,
								       const Sphere_point_3 &z) const{
  /*FT p[3];
  p[plane_coordinate(0).index()]= table_->center(x)[plane_coordinate(0).index()];
  p[plane_coordinate(1).index()]= table_->center(y)[plane_coordinate(1).index()];
  p[sweep_coordinate().index()]=0;
  FT v[3];
  v[plane_coordinate(0).index()]=0;
  v[plane_coordinate(1).index()]=0;
  v[sweep_coordinate().index()]=1;*/
  Line_3 l= z.line();//(Point_3(p[0], p[1], p[2]), Vector_3(v[0],v[1],v[2]));
    
  Sphere_point_3 p0(table_->sphere(s), l);

  if (!p0.is_valid()) {
    return CGAL::ON_UNBOUNDED_SIDE;
  } else {
    Sphere_point_3 p1(table_->sphere(s), 
		      l.opposite());
    CGAL::Comparison_result c0= compare_depths(p0, z);
    CGAL::Comparison_result c1= compare_depths(p1, z);
    if (c0 == CGAL::EQUAL || c1== CGAL::EQUAL) {
      //std::cout << "Equal" << std::endl;
      return CGAL::ON_BOUNDARY;
    }
    //CGAL_assertion(p0 <= p1);
    if (c0== CGAL::SMALLER && c1 == CGAL::LARGER){
      return CGAL::ON_BOUNDED_SIDE;
    } else {
      return CGAL::ON_UNBOUNDED_SIDE;
    }
  }
}



CGAL_AOS3_TEMPLATE
CGAL::Bounded_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::bounded_side_of_sphere(Sphere_3_key s,
								       Sphere_3_key x,
								       Sphere_3_key y,
								       const Sphere_point_3 &z) const {
  FT p[3];
  p[CGAL_AOS3_INTERNAL_NS::plane_coordinate(0).index()]= table_->center(x)[CGAL_AOS3_INTERNAL_NS::plane_coordinate(0).index()];
  p[CGAL_AOS3_INTERNAL_NS::plane_coordinate(1).index()]= table_->center(y)[CGAL_AOS3_INTERNAL_NS::plane_coordinate(1).index()];
  p[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()]=0;
  FT v[3];
  v[CGAL_AOS3_INTERNAL_NS::plane_coordinate(0).index()]=0;
  v[CGAL_AOS3_INTERNAL_NS::plane_coordinate(1).index()]=0;
  v[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()]=1;
  Line_3 l(Point_3(p[0], p[1], p[2]), Vector_3(v[0],v[1],v[2]));
    
  Sphere_point_3 p0(table_->sphere(s), l);

  if (!p0.is_valid()) {
    return CGAL::ON_UNBOUNDED_SIDE;
  } else {
    Sphere_point_3 p1(table_->sphere(s), 
		      l.opposite());
    CGAL::Comparison_result c0= compare_depths(p0, z);
    CGAL::Comparison_result c1= compare_depths(p1, z);
    if (c0 == CGAL::EQUAL || c1== CGAL::EQUAL) {
      //std::cout << "Equal" << std::endl;
      return CGAL::ON_BOUNDARY;
    }
    //CGAL_assertion(p0 <= p1);
    if (c0== CGAL::SMALLER && c1 == CGAL::LARGER){
      return CGAL::ON_BOUNDED_SIDE;
    } else {
      return CGAL::ON_UNBOUNDED_SIDE;
    }
  }
}


CGAL_AOS3_TEMPLATE
bool 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::is_over_circle_c(Sphere_3_key s,
								 const Sphere_point_3& z,
								 Coordinate_index C) const {
  Line_3 l= z.line();
  FT p[3], v[3];
  for (unsigned int i=0; i< 3; ++i) {
    p[i]= l.point()[i];
    v[i]= l.to_vector()[i];
  }
  Coordinate_index oc= CGAL_AOS3_INTERNAL_NS::other_plane_coordinate(C);
  p[oc.index()]= table_->center(s)[oc.index()];
  v[oc.index()]= 0;
  Line_3 nl=Line_3(Point_3(p[0], p[1], p[2]),Vector_3(v[0], v[1], v[2]));
  Sphere_point_3 sp(table_->sphere(s), nl);
  if (!sp.is_valid()) return false;
  else {
    Sphere_point_3 spo(table_->sphere(s), nl.opposite());
    CGAL::Comparison_result c= sp.compare(z, CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
    CGAL::Comparison_result co= spo.compare(z, CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
    if (co == c) return true;
    else return false;
  }
}

CGAL_AOS3_TEMPLATE
CGAL::Oriented_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::oriented_side_of_separating_plane(Sphere_3_key a,
										  Sphere_3_key b, 
										  const Sphere_point_3& ep) const {
  Plane_3 sp= table_->separating_plane(a,b);
  return oriented_side(sp, ep);
}


CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_depths(const Sphere_point_3 &a, 
							       const Sphere_point_3 &b) const{
  return a.compare(b, CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
}


CGAL_AOS3_TEMPLATE
CGAL::Oriented_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::oriented_side(const Plane_3 &p,
							      const Sphere_point_3 &s) const{
  CGAL::Object o= geom_traits_object().intersect_3_object()(p, s.line());
  Point_3 pt;
  Line_3 l;
  if (CGAL::assign(pt, o)){
    CGAL::Comparison_result c= s.compare_on_line(pt);
    CGAL::Sign sn= CGAL::sign(p.orthogonal_vector() * s.line().to_vector());
    if (sn ==CGAL::POSITIVE) return CGAL::enum_cast<CGAL::Oriented_side>(c);
    else return CGAL::enum_cast<CGAL::Oriented_side>(-c);
  } else if (CGAL::assign( l, o)){
    return CGAL::ON_ORIENTED_BOUNDARY;
  } else {
    return p.oriented_side(s.line().point());
  }
}

 
CGAL_END_NAMESPACE
