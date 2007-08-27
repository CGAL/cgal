#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Coordinate_index.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_key.h>
CGAL_BEGIN_NAMESPACE

#define CGAL_AOS3_CANONICAL_PREDICATE_LINE(l, sphere, t, equal_value, inside_value, below_value, above_value, miss_value) \
  std::cout << "Line is " << l << std::endl;				\
  Sphere_point_3 sp0(sphere, l);					\
  if (!sp0.is_valid()) {						\
      std::cout << "Line misses" << std::endl;				\
      miss_value;							\
  } else {								\
    Sphere_point_3 sp1(sphere, l.opposite());				\
    Coordinate_index sc=CGAL_AOS3_INTERNAL_NS::sweep_coordinate();	\
    Comparison_result c0= t.compare(sp0, sc);				\
    Comparison_result c1= t.compare(sp1, sc);				\
    if (c0 == CGAL::EQUAL || c1 == CGAL::EQUAL) {			\
      std::cout << "Equal" << std::endl;				\
      equal_value;							\
    } else if (c0 != c1) {						\
      std::cout << "Inside" << std::endl;				\
      inside_value;							\
    } else if (c0 == CGAL::SMALLER) {					\
      std::cout << "smaller" << std::endl;				\
      below_value;							\
    } else {								\
      std::cout << "larger" << std::endl;				\
      above_value;							\
    }									\
  }									\
  		
							
#define CGAL_AOS3_CANONICAL_PREDICATE(plane0, plane1, sphere, t, equal_value, inside_value, below_value, above_value, miss_value) \
  Plane_3 p0= plane0;							\
  std::cout << "Plane0 is " << p0 << std::endl;				\
  Plane_3 p1= plane1;							\
  std::cout << "Plane1 is " << p1 << std::endl;				\
  Object o= table_->geom_traits_object().intersect_3_object()(p0, p1);	\
  Line_3 l;								\
  if (assign(l,o) ){							\
    CGAL_AOS3_CANONICAL_PREDICATE_LINE(l, sphere, t, equal_value,	\
				       inside_value, below_value,	\
				       above_value, miss_value);	\
  } else {								\
    std::cout << "Planes miss" << std::endl;				\
    miss_value;								\
  }									

// predicates -------------------------------------------------------------------------------

CGAL_AOS3_TEMPLATE
bool Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::intersects(Sphere_3_key a,
								Sphere_3_key b) const {
  typedef CGAL_AOS3_TYPENAME Table::Equal_centers_exception SE;
  try {
    Point_3 eqpt= table_->equipower_point(a,b);
    std::cout << eqpt  << std::endl;
    CGAL::Bounded_side bs=table_->sphere(a).bounded_side(eqpt);
    //CGAL::Bounded_side obs=table_->sphere(b).bounded_side(eqpt);
    return bs != CGAL::ON_UNBOUNDED_SIDE;
  } catch (SE q) {
    std::cout << "Equal centers for " << a << " and " << b << std::endl;
    return false;
  }
  /*}
    return ret;*/
}



CGAL_AOS3_TEMPLATE
bool Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::intersects(Sphere_3_key a,
								Sphere_3_key b,
								Sphere_3_key c) const {
  CGAL_precondition(intersects(a,b));
  CGAL_precondition(intersects(b,c));
  CGAL_precondition(intersects(c,a));
  // NOTE could be better, I think
  Plane_3 pab = table_->equipower_plane(a,b);
  Plane_3 pac = table_->equipower_plane(a,c);
  CGAL::Object o= geom_traits_object().intersect_3_object()(pab, pac);
  Line_3 l;
  if (CGAL::assign(l,o)){
    Sphere_point_3 sp(table_->sphere(a), l);
    return sp.is_valid();
  } else {
    std::cout << "Degeneracy for " << a << " " << b << " " << c
	      << " with planes " << pac << " and " << pab << std::endl;
    return pab== pac || pab.opposite() == pac;
  }
}


#if 0

CGAL_AOS3_TEMPLATE
bool 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::center_is_in_slab_c(Sphere_3_key sphere_for_slab, 
								    const Event_point_3 &point,
								    const Sphere_3_key sphere_for_point,
								    Coordinate_index C) const {
  FT dc2= CGAL::square(table_->center(sphere_for_slab)[C.index()]
		       - table_->center(sphere_for_point)[C.index()]);
  FT R2= table_->squared_radius(sphere_for_slab);
  if (dc2 > R2) return false;
  FT dt2= CGAL::square(table_->center(sphere_for_slab)[CGAL_AOS3_INTERNAL_NS::Sweep_coordinate::index()]
		       - table_->center(sphere_for_point)[CGAL_AOS3_INTERNAL_NS::Sweep_coordinate::index()]);
  
  FT r2= R2-dt2;
  bool myret= (r2 > dc2);
  CGAL_assertion(myret== point_is_in_slab_c(sphere_for_slab,
					    sphere_events(sphere_for_point).first,
					    C));
  return myret;
}


/*
  Project line onto center plane for sphere
  if it misses sphere, then we are outside plane

  compare plane intersection point (i.e. input) to slab sphere intersection point
*/
CGAL_AOS3_TEMPLATE
Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_point_to_circle_c(const Event_point_3 &point,
									  Sphere_3_key sphere_for_slab, 
									  Coordinate_index C) const {
  FT lpt[3]={point.line().point().x(), point.line().point().y(), point.line().point().z()};
  lpt[other_plane_coordinate(C).index()]=table_->center(sphere_for_slab)[other_plane_coordinate(C).index()];
  FT lv[3]={point.line().to_vector().x(), point.line().to_vector().y(), point.line().to_vector().z()};
  lv[other_plane_coordinate(C).index()]=0;
  Line_3 l3p(Point_3(lpt[0], lpt[1], lpt[2]),
	     Vector_3(lv[0], lv[1], lv[2]));
  Event_point_3 nsp(table_->sphere(sphere_for_slab), l3p);
  if (nsp.is_valid()) {
    Event_point_3 onsp(table_->sphere(sphere_for_slab), l3p.opposite());
    if (nsp > onsp) std::swap(nsp, onsp);
    return (point >= nsp && point <= onsp);
  } else {
    return false;
  }
  
  // if point is above center and intersection then it is outside? 
  
}
#endif



CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_point_to_circle_circle_c(const Sphere_point_3 &pt,
										 Sphere_3_key sphere0,
										 Sphere_3_key sphere1,
										 Coordinate_index C) const{
 


  Plane_3 eqp=table_->equipower_plane(sphere0, sphere1);
  CGAL::Comparison_result cr= compare_point_to_equipower_line_c(pt, sphere0, sphere1, C);
  //Oriented_side ossp= center_oriented_side_of_separating_plane_c(k, sphere0, sphere1, C);
  Plane_3 sp= table_->separating_plane(sphere0, sphere1);
  CGAL::Comparison_result epqc= compare(sp.orthogonal_vector()[C.index()],0);
  if (cr != epqc ) {
    return cr;
  }
  
  CGAL_AOS3_CANONICAL_PREDICATE(eqp,
				const_c_plane(pt, C),
				table_->sphere(sphere0),
				pt,
				return CGAL::EQUAL,
				return Comparison_result(-cr),
				return cr,
				return cr,
				return cr);
}

CGAL_AOS3_TEMPLATE
CGAL::Comparison_result
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_center_to_circle_circle_c(Sphere_3_key k,
										  const Event_point_3 &t,
										  Sphere_3_key sphere0,
										  Sphere_3_key sphere1,
										  Coordinate_index C) const{
  //std::cout << "EQP is " << table_->equipower_point(sphere0, sphere1) << " 1/10000 1" << std::endl;

  //std::cout << "EQPlane is " << table_->equipower_plane(sphere0, sphere1)<< std::endl;
  //std::cout << "Sep plane is " << table_->separating_plane(sphere0, sphere1)<< std::endl;
 
  CGAL::Comparison_result cr= compare_center_to_equipower_line_c(k, t, sphere0, sphere1, C);
  //Oriented_side ossp= center_oriented_side_of_separating_plane_c(k, sphere0, sphere1, C);
  Plane_3 sp= table_->separating_plane(sphere0, sphere1);
  CGAL::Comparison_result epqc= compare(sp.orthogonal_vector()[C.index()],0);
  std::cout << "comps are " << cr << " " << epqc << std::endl;
  if (cr != epqc ) {
    return cr;
  }
  Plane_3 eqp=table_->equipower_plane(sphere0, sphere1);

 
  CGAL_AOS3_CANONICAL_PREDICATE(eqp,
				rule_plane(k, C),
				table_->sphere(sphere0),
				t,
				return CGAL::EQUAL,
				return Comparison_result(-cr),
				return cr,
				return cr, 
				return cr);
}


CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_point_to_circle_rule_c(const Sphere_point_3 &pt,
									       Sphere_3_key sphere,
									       Sphere_3_key rule,
									       Coordinate_index rule_coordinate,
									       Coordinate_index C) const{
  if (C== rule_coordinate) return compare_point_to_rule_c(pt, rule, rule_coordinate);

  Comparison_result crrs= compare_sphere_centers_c(rule, sphere, C);
  Comparison_result crcs= compare_point_to_rule_c(pt, sphere, C);
  CGAL_assertion(crrs != CGAL::EQUAL);
  if (crrs != crcs) return crcs;

  //Comparison_result cr= compare_point_to_rule_c(pt, rule, C);

  //Comparison_result cr= compare_point_to_rule_c(pt, sphere, C);
  CGAL_AOS3_CANONICAL_PREDICATE(const_c_plane(pt, C),
				rule_plane(rule, rule_coordinate),
				table_->sphere(sphere),
				pt,
				return CGAL::EQUAL,
				return Comparison_result(-crrs),
				return crrs,
				return crrs,
				return crrs);
}
  
CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_center_to_circle_rule_c(Sphere_3_key center,
										const Sphere_point_3 &t,
										Sphere_3_key sphere,
										Sphere_3_key rule,
										Coordinate_index rule_coordinate,
										Coordinate_index C) const{
  if (C== rule_coordinate) return CGAL::compare(table_->center(center)[rule_coordinate.index()],
						table_->center(rule)[rule_coordinate.index()]);
  Comparison_result crrs= compare_sphere_centers_c(rule, sphere, C);
  Comparison_result crcs= compare_sphere_centers_c(center, sphere, C);
  CGAL_assertion(crrs != CGAL::EQUAL);
  if (crrs != crcs) return crcs;

  //Comparison_result cr= compare_sphere_centers_c(center, rule, C);
  //if (cr == crrs) return cr;

  CGAL_AOS3_CANONICAL_PREDICATE(rule_plane(center, C),
				rule_plane(rule, rule_coordinate),
				table_->sphere(sphere), 
				t,
				return CGAL::EQUAL,
					return Comparison_result(-crrs),
				return crrs,
				return crrs,
				return crrs);
}
  
CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_point_to_rule_c(const Sphere_point_3 &pt,
									Sphere_3_key rule,
									Coordinate_index C) const{
  return pt.compare(table_->center(rule)[C.index()], C);
}





CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_point_to_equipower_line_c(const Sphere_point_3 &sp,
										   Sphere_3_key a, 
										   Sphere_3_key b,
										   Coordinate_index C) const{
  // NOTE redo this with computing it directly
  CGAL_assertion(C!=CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
  Line_3 l= table_->equipower_line(a,b);
  const Line_3 &plt= sp.line();
  Line_2 l2(Point_2(l.point()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()],
		    l.point()[C.index()]),
	    Vector_2(l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()],
		     l.to_vector()[C.index()]));
  Line_2 plt2(Point_2(plt.point()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()],
		      plt.point()[C.index()]),
	      Vector_2(plt.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()],
		       plt.to_vector()[C.index()]));
  //std::cout << "l2 is " << l2 << std::endl;
  //std::cout << "pl2 is " << plt2 << std::endl;
  Object o= geom_traits_object().intersect_2_object()(l2, plt2);
  Point_2 p;
  FT cmp=0;
  Comparison_result c0= compare(plt2.y_at_x(cmp), l2.y_at_x(cmp));
  if (c0 == EQUAL) {
    cmp=1;
    c0=  compare(plt2.y_at_x(cmp), l2.y_at_x(cmp));
    if (c0 == EQUAL) {
      return EQUAL;
    }
  }
  if (assign(p,o)) {

    //std::cout << "Sign changes at " << p.x() << std::endl;
    if (sp.compare(p.x(),CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) == EQUAL) return EQUAL;
    if (p.x()>0 && sp.compare(p.x(), CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) == LARGER
	|| p.x()<0 && sp.compare(p.x(), CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) == SMALLER ){
      return Comparison_result(-c0);
    } else return c0;
  } else {
    return c0;
  }
}

CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_center_to_equipower_line_c(Sphere_3_key k,
										   const Event_point_3 &t,
										   Sphere_3_key a, 
										   Sphere_3_key b,
										   Coordinate_index C) const{
  // NOTE redo this with computing it directly
    // NOTE redo this with computing it directly
  CGAL_assertion(C!=CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
  Line_3 l= table_->equipower_line(a,b);
  //std::cout << "EQPlane is " << table_->equipower_plane(a,b)<< std::endl;
  //std::cout << "EQLineis " << table_->equipower_line(a,b)<< std::endl;
   
  Line_2 l2(Point_2(l.point()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()],
		    l.point()[C.index()]),
	    Vector_2(l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()],
		     l.to_vector()[C.index()]));
  Line_2 plt2(Point_2(0,
		      table_->center(k)[C.index()]),
	      Vector_2(1,0));
  //std::cout << "l2 is " << l2 << std::endl;
  //std::cout << "pl2 is " << plt2 << std::endl;
  CGAL_assertion(plt2.y_at_x(0) == table_->center(k)[C.index()]);

  Object o= geom_traits_object().intersect_2_object()(l2, plt2);
  Point_2 p;
  FT cmp=0;
  Comparison_result c0= compare(table_->center(k)[C.index()], l2.y_at_x(cmp));
  if (c0 == EQUAL) {
    cmp=1;
    //std::cout << "falling back" << std::endl;
    c0=  compare(table_->center(k)[C.index()], l2.y_at_x(cmp));
    if (c0 == EQUAL) {
      return EQUAL;
    }
  }
  if (assign(p,o)) {

    //std::cout << "Sign changes at " << p.x() << std::endl;
    if (t.compare(p.x(),CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) == EQUAL) return EQUAL;
    if (p.x()>cmp && t.compare(p.x(), CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) == LARGER
	|| p.x()<cmp  && t.compare(p.x(), CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) == SMALLER ){
      //std::cout << " between" << std::endl;
      return Comparison_result(-c0);
    } else return c0;
  } else {
    return c0;
  }
}




/*CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_point_to_equipower_point_c(const Sphere_point_3 &sp,
										   Sphere_3_key a, 
										   Sphere_3_key b,
										   Coordinate_index C) const{
  // NOTE redo this with computing it directly
  Point_3 pt= table_->equipower_point(a,b);
  return sp.compare(pt, C);
}

CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_center_to_equipower_point_c(Sphere_3_key k,
										   Sphere_3_key a, 
										   Sphere_3_key b,
										   Coordinate_index C) const{
  // NOTE redo this with computing it directly
  Point_3 pt= table_->equipower_point(a,b);
  return compare(table_->center(k)[C.index()], pt[C.index()]);
}
*/



CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_points_c(const Sphere_point_3 &a,
								 const Sphere_point_3 &b,
								 Coordinate_index C) const{
  return a.compare(b, C);
}
										   

CGAL_AOS3_TEMPLATE
CGAL::Bounded_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::point_bounded_side_of_sphere(const Sphere_point_3 &pt,
									     Sphere_3_key s
									     ) const{
  
  // This is a bit silly but makes it uniform
  CGAL_AOS3_CANONICAL_PREDICATE(const_c_plane(pt, CGAL_AOS3_INTERNAL_NS::plane_coordinate(0)),
				const_c_plane(pt, CGAL_AOS3_INTERNAL_NS::plane_coordinate(1)),
				table_->sphere(s),
				pt,
				CGAL_assertion(l== pt.line() || l.opposite() ==pt.line()); return ON_BOUNDARY,
				CGAL_assertion(l== pt.line() || l.opposite() ==pt.line()); return ON_BOUNDED_SIDE,
				CGAL_assertion(l== pt.line() || l.opposite() ==pt.line()); return ON_UNBOUNDED_SIDE,
				CGAL_assertion(l== pt.line() || l.opposite() ==pt.line()); return ON_UNBOUNDED_SIDE,
				CGAL_assertion(l== pt.line() || l.opposite() ==pt.line()); return ON_UNBOUNDED_SIDE);
  /*
  const Line_3 &l= pt.line();//(Point_3(p[0], p[1], p[2]), Vector_3(v[0],v[1],v[2]));
    
  Sphere_point_3 p0(table_->sphere(s), l);

  if (!p0.is_valid()) {
    return CGAL::ON_UNBOUNDED_SIDE;
  } else {
    Sphere_point_3 p1(table_->sphere(s), 
		      l.opposite());
    CGAL::Comparison_result c0= pt.compare(p0, CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
    CGAL::Comparison_result c1= pt.compare(p1, CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
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
    }*/
}

CGAL_AOS3_TEMPLATE
CGAL::Bounded_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::center_bounded_side_of_sphere(const Event_point_3 &t,
									      Sphere_3_key k,
									      Sphere_3_key s
									      ) const{
  
  // This is a bit silly but makes it uniform
  CGAL_AOS3_CANONICAL_PREDICATE(rule_plane(k, CGAL_AOS3_INTERNAL_NS::plane_coordinate(0)),
				rule_plane(k, CGAL_AOS3_INTERNAL_NS::plane_coordinate(1)),
				table_->sphere(s),
				t,
				return ON_BOUNDARY,
				return ON_BOUNDED_SIDE,
				return ON_UNBOUNDED_SIDE,
				return ON_UNBOUNDED_SIDE,
				return ON_UNBOUNDED_SIDE);
  /*
  const Line_3 &l= pt.line();//(Point_3(p[0], p[1], p[2]), Vector_3(v[0],v[1],v[2]));
    
  Sphere_point_3 p0(table_->sphere(s), l);

  if (!p0.is_valid()) {
    return CGAL::ON_UNBOUNDED_SIDE;
  } else {
    Sphere_point_3 p1(table_->sphere(s), 
		      l.opposite());
    CGAL::Comparison_result c0= pt.compare(p0, CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
    CGAL::Comparison_result c1= pt.compare(p1, CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
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
    }*/
}


CGAL_AOS3_TEMPLATE
CGAL::Bounded_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::center_bounded_side_of_circle_c(Sphere_3_key k,
										const Event_point_3 &t,
										Sphere_3_key s,
										Coordinate_index C) const{
  //Comparison_result cr= compare_sphere_centers_c(k,s, C);
  CGAL_AOS3_CANONICAL_PREDICATE(rule_plane(s, other_plane_coordinate(C)),
				rule_plane(k, C),
				table_->sphere(s),
				t,
				return ON_BOUNDARY,
				return ON_BOUNDED_SIDE,
				return ON_UNBOUNDED_SIDE,
				return ON_UNBOUNDED_SIDE,
				return ON_UNBOUNDED_SIDE);
}


CGAL_AOS3_TEMPLATE
CGAL::Bounded_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::point_bounded_side_of_sphere_c(const Sphere_point_3 &pt,
									       Sphere_3_key s,
									       Coordinate_index C) const{
 
  if (C== CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) {
    // also silly but uniform

    CGAL_AOS3_CANONICAL_PREDICATE(rule_plane(s, CGAL_AOS3_INTERNAL_NS::plane_coordinate(0)),
				  rule_plane(s, CGAL_AOS3_INTERNAL_NS::plane_coordinate(1)),
				  table_->sphere(s),
				  pt,
				  return ON_BOUNDARY,
				  return ON_BOUNDED_SIDE,
				  return ON_UNBOUNDED_SIDE,
				  return ON_UNBOUNDED_SIDE,
				  return ON_UNBOUNDED_SIDE);
    /*Event_pair ep= sphere_events(s);
    if (pt.compare(ep.first, C)==EQUAL || pt.compare(ep.second, C)==EQUAL) return CGAL::ON_BOUNDARY;
    else if (pt.compare(ep.first, C) == LARGER && pt.compare(ep.second,C)==SMALLER) return CGAL::ON_BOUNDED_SIDE;
    else return CGAL::ON_UNBOUNDED_SIDE;*/
  } else {
    CGAL_AOS3_CANONICAL_PREDICATE(const_c_plane(pt, C),
				  rule_plane(s, CGAL_AOS3_INTERNAL_NS::other_plane_coordinate(C)),
				  table_->sphere(s),
				  pt,
				  return ON_BOUNDARY,
				  return ON_BOUNDED_SIDE,
				  return ON_UNBOUNDED_SIDE,
				  return ON_UNBOUNDED_SIDE,
				  return ON_UNBOUNDED_SIDE);
  }
}

CGAL_AOS3_TEMPLATE
CGAL::Bounded_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::rules_bounded_side_of_sphere(const Sphere_point_3 &t,
									     Sphere_3_key x,
									     Sphere_3_key y,
									     Sphere_3_key sphere) const {
  CGAL_AOS3_CANONICAL_PREDICATE(rule_plane(x, CGAL_AOS3_INTERNAL_NS::plane_coordinate(0)),
				rule_plane(y, CGAL_AOS3_INTERNAL_NS::plane_coordinate(1)),
				table_->sphere(sphere),
				t,
				return ON_BOUNDARY,
				return ON_BOUNDED_SIDE,
				return ON_UNBOUNDED_SIDE,
				return ON_UNBOUNDED_SIDE,
				return ON_UNBOUNDED_SIDE);
}


CGAL_AOS3_TEMPLATE
CGAL::Comparison_result
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_sphere_centers_c(Sphere_3_key a, Sphere_3_key b,
									 Coordinate_index C) const {
  return CGAL::compare(table_->center(a)[C.index()], table_->center(b)[C.index()]);
}



CGAL_AOS3_TEMPLATE
CGAL::Oriented_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::point_oriented_side_of_separating_plane( const Sphere_point_3& ep,
											 Sphere_3_key a,
											 Sphere_3_key b ) const {
  Plane_3 sp= table_->separating_plane(a,b);
  return oriented_side(sp, ep);
}

CGAL_AOS3_TEMPLATE
CGAL::Oriented_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::center_oriented_side_of_separating_plane(Sphere_3_key center,
											 const Sphere_point_3& ep,
											 Sphere_3_key a,
											 Sphere_3_key b ) const {
  CGAL::Object o= geom_traits_object().intersect_3_object()(table_->separating_plane(a,b),
							    Line_3(table_->center(center),
								   CGAL_AOS3_INTERNAL_NS::sweep_vector<Vector_3>()));
  Oriented_side os= table_->separating_plane(a,b).oriented_side(table_->center(center));
  Point_3 pt;
  Line_3 l;
  if (CGAL::assign(pt, o)){
    Comparison_result cc= ep.compare(table_->center(center),CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
    Comparison_result ic= ep.compare(pt,CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
    if (ic == CGAL::EQUAL) return CGAL::ON_ORIENTED_BOUNDARY;
    else if (cc == ic) return os;
    else return Oriented_side(-os);
  } else if (CGAL::assign( l, o)){
    return CGAL::ON_ORIENTED_BOUNDARY;
  } else {
    return os;
  }
}




// Internal predicates




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

























#if 0



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
CGAL::Bounded_side 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG
::bounded_side_of_sphere_projected(const Sphere_point_3 &t,
				   Sphere_3_key sphere,
				   Sphere_3_key oplane,
				   const Sphere_point_3 &pt,
				   Coordinate_index C) const{
  FT rd[3]={0,0,0};
  rd[other_plane_coordinate(C).index()]=1;
  Plane_3 p(pt.line().point(),
	    pt.line().point() + pt.line().to_vector(),
	    pt.line().point() + Vector_3(rd[0], rd[1], rd[2]));
  std::cout << "Rule plane is " << p << std::endl;
  FT v[3]={0,0,0};
  v[other_plane_coordinate(C).index()]=1;
  Plane_3 rp(table_->center(oplane), Vector_3(v[0], v[1], v[2]));
  std::cout << "Other rule plane is " << rp << std::endl;
  Line_3 l;
  CGAL::Object o= geom_traits_object().intersect_3_object()(rp, p);
  if (!CGAL::assign(l,o)) {
    CGAL_assertion(0);
    return CGAL::ON_UNBOUNDED_SIDE;
  } else {
    std::cout << "Line is " << l << std::endl;
    Sphere_point_3 p0(table_->sphere(sphere), l);
    if (!p0.is_valid()) {
      std::cout << "Line misses" << std::endl;
      return CGAL::ON_UNBOUNDED_SIDE;
    } else {
      Sphere_point_3 p1(table_->sphere(sphere), l.opposite());
      std::cout << "Points are " << p0 << " " << p1 << std::endl;
      CGAL::Comparison_result c0= compare_depths(t, p0), c1= compare_depths(t, p1);
      std::cout << "Depths are " << c0 << " " << c1 << std::endl;
      if (c0==CGAL::EQUAL || c1 == CGAL::EQUAL) {
	return CGAL::ON_BOUNDARY;
      } else if (c0 != c1) {
	return CGAL::ON_BOUNDED_SIDE;
      } else {
	return CGAL::ON_UNBOUNDED_SIDE;
      }
    }
  }
}







CGAL_AOS3_TEMPLATE
bool Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::sphere_intersects_rule(Sphere_3_key s, 
									    Sphere_3_key rule, 
									    Coordinate_index C) const {
  FT d= table_->center(s)[C.index()]-table_->center(rule)[C.index()];
  if (CGAL::square(d) <= table_->sphere(s).squared_radius()) return true;
  else return false;
}

CGAL_AOS3_TEMPLATE
bool Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::sphere_over_rule(Sphere_3_key s, 
								      const Sphere_point_3 &sp,
								      Coordinate_index C) const {
  FT v[3]={0,0,0};
  CGAL::Comparison_result cr=sp.compare(table_->center(s), C);
  if (cr == CGAL::SMALLER) {
    v[C.index()]=-1;
  } else if (cr == CGAL::LARGER) {
    v[C.index()]=1;
  } else {
    return true;
  }
  Sphere_point_3 ep(table_->sphere(s), Line_3(table_->center(s), 
					      Vector_3(v[0], v[1], v[2])));
  return ep.compare(sp, C) != cr;
}




CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::compare_sphere_sphere_at_sweep(const Sphere_point_3 &t,
									       Sphere_3_key sphere0,
									       Sphere_3_key sphere1,
									       const Sphere_point_3 &pt,
									       Coordinate_index C) const {
  {
    std::cout << "Compare SS at sweep " << t << " spheres " 
	      <<  sphere0 << " " << sphere1 
	      << " point " << pt << " line " << pt.line() << " coord " << C 
	      << std::endl;

    // total equality is actually quite easy
    /*Sphere_point_3 pt0(table_->sphere(sphere0), pt.line());
      if (pt0.compare(pt, C) == CGAL::EQUAL && ) {
      Sphere_point_3 pt1(table_->sphere(sphere1), pt.line());
      if (pt1.compare(pt, C) == CGAL::EQUAL) return CGAL::EQUAL;
      }*/

    // compare to the extrema of the two spheres to makes sure that the point overlaps them
    bool i0= sphere_over_rule(sphere0, pt, C);
     
     
    if (!i0) {
      CGAL::Comparison_result sc0= compare_sphere_center_c(sphere0, pt, C);
      std::cout << "Returning " << sc0 << " since sphere 0 misses" << std::endl;
      return sc0;
    }

    bool i1= sphere_over_rule(sphere1, pt, C);

    if (!i1) {
      CGAL::Comparison_result sc1= compare_sphere_center_c(sphere1, pt, C);
      std::cout << "Returning " << sc1 << " since sphere 0 misses" << std::endl;
      return sc1;
    }
  }


  Plane_3 eqp= table_->equipower_plane(sphere0, sphere1);
  std::cout << "EQP is " << eqp << std::endl;
   
  FT v[3]={0,0,0};
  v[CGAL_AOS3_INTERNAL_NS::other_plane_coordinate(C).index()]=1;
   
  Plane_3 rule_plane(pt.line().point(),
		     pt.line().point()+ Vector_3(v[0], v[1], v[2]),
		     pt.line().point()+ pt.line().to_vector());
  if (rule_plane.orthogonal_vector()[C.index()] < 0) {
    rule_plane=rule_plane.opposite();
  }
  std::cout << "rule_plane is " << rule_plane << std::endl;
  CGAL_assertion(rule_plane != Plane_3());
  Line_3 l;
  CGAL::Object o= geom_traits_object().intersect_3_object()(eqp, rule_plane);
  if (!CGAL::assign(l,o)) {
    CGAL::Comparison_result cr= compare_equipower_point_c(sphere0, sphere1,
							  pt, C);
    std::cout << "Returning " << cr << " since sphere the planes miss" << std::endl;
    return cr;
  }
 
  std::cout << "The intersection line is " << l << std::endl;
  

   
  Sphere_point_3 switch_point_0(table_->sphere(sphere0), l);
   
  if (!switch_point_0.is_valid()){
     
    CGAL::Comparison_result cr=  compare_equipower_point_c(sphere0,sphere1,
							   pt, C);
    std::cout << "Returning " << cr << " line misses" << std::endl;
    return cr;
  }
 
  //if (switch_point_0 > switch_point_1) {
  //  std::swap(switch_point_0, switch_point_1);
  //}
  Plane_3 sepp= table_->separating_plane(sphere0, sphere1);
  std::cout << "separating_plane is " << sepp << std::endl;
   
  Sphere_point_3 frontp= intersection_2_events(sphere0, sphere1).first;
  std::cout << "The front point is " << frontp << std::endl;
  CGAL::Oriented_side fpos= oriented_side(rule_plane, frontp);
  bool flip=false;
  if (fpos == CGAL::ON_ORIENTED_BOUNDARY) {
    std::cout << "Front point is on boundary using normals ";
    fpos = CGAL::Oriented_side(CGAL::sign(sepp.orthogonal_vector()
					  * rule_plane.orthogonal_vector()));
    std::cout << "got " << fpos << std::endl;
  } 

  std::cout << "Oriented_side is " << fpos << std::endl;
  CGAL::Comparison_result spcr0= switch_point_0.compare(t, 
							CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
   
  Sphere_point_3 switch_point_1(table_->sphere(sphere0), l.opposite());
  std::cout << "switch points are " << switch_point_0 
	    << " and " << switch_point_1 << std::endl;
  CGAL::Comparison_result spcr1= switch_point_1.compare(t, 
							CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
  std::cout << "spcr is " << spcr0 << spcr1 << std::endl;
   
   
  CGAL::Oriented_side os_0= oriented_side(sepp, switch_point_0);
  CGAL::Oriented_side os_1= oriented_side(sepp, switch_point_1);
  std::cout << "os are " << os_0 << os_1 << std::endl;


  if (spcr0 == CGAL::EQUAL && os_0 == CGAL::ON_POSITIVE_SIDE
      || spcr1 == CGAL::EQUAL && os_1 == CGAL::ON_POSITIVE_SIDE) {
    std::cout << "Returning equal" << std::endl;
    return CGAL::EQUAL;
  }
   
  if (os_0 == CGAL::ON_POSITIVE_SIDE && os_1 == CGAL::ON_POSITIVE_SIDE) {
    flip = (spcr0 != spcr1);
  } else if (os_0 == CGAL::ON_POSITIVE_SIDE) {
    flip = (spcr0 == CGAL::SMALLER);
  } else if (os_1 == CGAL::ON_POSITIVE_SIDE) {
    flip = (spcr1 == CGAL::LARGER);
  }

  if (!flip) {
    return CGAL::Comparison_result(fpos);
  } else {
    return CGAL::Comparison_result(-fpos);
  }
}





#endif


// Constructions----------------------------------------------------------------------------------------






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
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_pair 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::intersection_2_events(Sphere_3_key a, Sphere_3_key b) const {
  std::cout << "Intersection events for " << a << " and " << b << std::endl;
  if (!intersects(a,b)) {
    return Event_pair();
  }
  Plane_3 eqp= table_->equipower_plane(a,b);
  Plane_3 sp= table_->separating_plane(a,b);
  //CGAL_assertion(eqp.has_on(eqpoint));
  CGAL::Object o= geom_traits_object().intersect_3_object()(eqp, sp);
  std::cout << "eqp= " << eqp << std::endl << "sp= " << sp << std::endl;
  //<< "eqpoint= " << eqpoint << std::endl;
  Line_3 l;
  if (!CGAL::assign(l,o)){
    std::cout << "Two spheres with same center" << std::endl;
    Event_point_3 ep(table_->equipower_point(a,b));
    return Event_pair(ep, Event_point_3());
  }
  Line_3 lf, lb;
  
  if (l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()] > 0){
    lf=l;
    lb=l.opposite();
  } else {
    lf=l.opposite();
    lb=l;
  }
  if (l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()]==0) {
    std::cerr << "Degeneracy, picking one for start of circle " 
	      << a << " " << b << std::endl;
  }
  
  std::cout << table_->sphere(a) << std::endl 
	    << "lf= " << lf << std::endl 
	    << "lb= " << lb << std::endl;
  Sphere_point_3 spa(table_->sphere(a), lf);
  Sphere_point_3 spb(table_->sphere(a), lb);
  std::cout << spa << std::endl << spb << std::endl;
  Event_pair ret(spa, spb);
  //CGAL_assertion(compare_depths(ret.first, ret.second) != CGAL::LARGER);
  return ret;
}
 












CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Sphere_3_key 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::new_sphere_3(const Sphere_3 &s) {
  return table_->new_sphere(s);
}




CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_pair 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::intersection_3_events(Sphere_3_key a,
								      Sphere_3_key b, 
								      Sphere_3_key c) const {
  CGAL_precondition(a != b && a != c && c!= a);
  if (!intersects(a,b, c)) {
    return Event_pair();
  }
  Plane_3 eqpab= table_->equipower_plane(a,b);
  Plane_3 eqpac= table_->equipower_plane(a,c);

  CGAL::Object o= geom_traits_object().intersect_3_object()(eqpab, eqpac);
  Line_3 l;
  if (!CGAL::assign(l,o)){
    // they intersect on the same circle. Ick.
    CGAL_assertion(0);
  }
  Line_3 lf, lb;
  CGAL_assertion(l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()] != 0);
  if (l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()] > 0){
    lf=l;
    lb=l.opposite();
  } else {
    lf=l.opposite();
    lb=l;
  }
  if (l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()] ==0) {
    std::cerr << "Degeneracy, picking one for start of threesome " 
	      << a << " " << b << " " << c << std::endl;
  }
  
  Event_pair ret(Sphere_point_3(table_->sphere(a), lf),
		 Sphere_point_3(table_->sphere(a), lb));
  //CGAL_assertion(compare_depths(ret.first, ret.second) != CGAL::LARGER);
  return ret;
}



CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Plane_3 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::rule_plane(Sphere_3_key r, Coordinate_index C) const {
  FT n[3]={0,0,0};
  n[C.index()]= 1;
  return Plane_3(table_->center(r), Vector_3(n[0], n[1], n[2]));
}




CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Plane_3 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::const_c_plane(const Sphere_point_3 &pt, Coordinate_index C) const {
  FT v[3]={0,0,0};
  v[other_plane_coordinate(C).index()]=1;

 return Plane_3(pt.line().point(), pt.line().point() + pt.line().to_vector(),
		pt.line().point()+ Vector_3(v[0], v[1], v[2]));
}



CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_pair 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::sphere_intersect_extremum_events(Sphere_3_key a, 
										 Coordinate_index C,
										 Sphere_3_key b) const {
  CGAL_precondition(intersects(a,b));
  Plane_3 eqpab= table_->equipower_plane(a,b);
  Plane_3 rp= rule_plane(a,C);
  CGAL::Object o= geom_traits_object().intersect_3_object()(eqpab, rp);
  Line_3 l;
  if (!CGAL::assign(l,o)){
    return Event_pair();
  } else {
    Sphere_point_3 sa(table_->sphere(a), l);
    if (!sa.is_valid()) return Event_pair();
    Event_point_3 pa(sa);
    Event_point_3 pb(Sphere_point_3(table_->sphere(a), l.opposite()));
    if (pb < pa) std::swap(pa, pb);
    return Event_pair(pa, pb);
  }
}


CGAL_AOS3_TEMPLATE
void 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::advance_circle_cross_rule_event(Sphere_3_key a, 
										Sphere_3_key b,
										Sphere_3_key rs,
										Coordinate_index C) {
  CGAL_AOS3_INTERNAL_NS::Sphere_key_upair k(a,b);
  CGAL_AOS3_TYPENAME Table::Event_pair_data &ed
    = table_->upair_data(k).cxr_events_[CGAL_AOS3_TYPENAME Table::UPair_data::KC_pair(rs, C)];
  ++ed.index_;
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_point_3 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::circle_cross_rule_event(Sphere_3_key a, 
									Sphere_3_key b,
									Sphere_3_key rs,
									Coordinate_index C) const {
  CGAL_AOS3_INTERNAL_NS::Sphere_key_upair k(a,b);
  CGAL_AOS3_TYPENAME Table::Event_pair_data &ed
    = table_->upair_data(k).cxr_events_[CGAL_AOS3_TYPENAME Table::UPair_data::KC_pair(rs, C)];

  if (ed.index_==-1) {
    Event_pair ep= circle_cross_rule_events(a,b,rs,C);
    ed.events_[0]=ep.first;
    ed.events_[1]=ep.second;
    ed.index_=0;
  }
  if (ed.index_ < 2) return ed.events_[ed.index_];
  else return Event_point_3();
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_pair 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::circle_cross_rule_events(Sphere_3_key a, 
									 Sphere_3_key b,
									 Sphere_3_key rs,
									 Coordinate_index C) const {
  Plane_3 eqp= table_->equipower_plane(a,b);
  /*FT na[3]={0,0,0};
    na[C.index()]=1;*/
  Plane_3 rp = rule_plane(rs, C); //(table_->center(rs), Vector_3(na[0], na[1], na[2]));
  CGAL::Object o=geom_traits_object().intersect_3_object()(eqp, rp);
  Line_3 l;
  if (CGAL::assign(l,o)) {
    Sphere_point_3 sp(table_->sphere(a), l);
    //Event_point_3 sp();
    if (!sp.is_valid()) {
      return Event_pair();
    } else {
      Event_point_3 nep(sp);
      Event_point_3 spp(Sphere_point_3(table_->sphere(a), l.opposite()));
      Event_pair ep(std::min(nep, spp), std::max(nep, spp));
      return ep;
    }
  } else {
    return Event_pair();
  }
}

CGAL_AOS3_TEMPLATE
void 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::advance_sphere_intersect_rule_rule_event(Sphere_3_key s, 
											 Sphere_3_key x,
											 Sphere_3_key y) {
  CGAL_AOS3_INTERNAL_NS::Sphere_key_triple k(s,x,y);
  CGAL_AOS3_TYPENAME Table::Event_pair_data &ed= table_->triple_data(k).srr_events_;
  ++ed.index_;
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_point_3 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::sphere_intersect_rule_rule_event(Sphere_3_key s, 
										 Sphere_3_key x,
										 Sphere_3_key y) const {
  CGAL_AOS3_INTERNAL_NS::Sphere_key_triple k(s,x,y);
  CGAL_AOS3_TYPENAME Table::Event_pair_data &ed= table_->triple_data(k).srr_events_;
  if (ed.index_==-1) {
    Event_pair ep= sphere_intersect_rule_rule_events(s,x,y);
    ed.events_[0]=ep.first;
    ed.events_[1]=ep.second;
    ed.index_=0;
  }
  if (ed.index_ < 2) return ed.events_[ed.index_];
  else return Event_point_3();
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::Event_pair 
Arrangement_of_spheres_traits_3 CGAL_AOS3_TARG::sphere_intersect_rule_rule_events(Sphere_3_key s, 
										  Sphere_3_key x,
										  Sphere_3_key y) const {
  CGAL_precondition(s != Sphere_3_key());
  CGAL_precondition(x != Sphere_3_key());
  CGAL_precondition(y != Sphere_3_key());

  std::cout << "SRR event for " << s << " " << x << " " << y << std::endl;

  Plane_3 x_plane= rule_plane(x, CGAL_AOS3_INTERNAL_NS::plane_coordinate(0));
  Plane_3 y_plane= rule_plane(y, CGAL_AOS3_INTERNAL_NS::plane_coordinate(1));
  CGAL::Object o= geom_traits_object().intersect_3_object()(x_plane, y_plane);
  Line_3 l;
  if (!CGAL::assign(l,o)){
    return Event_pair();
  } else {
    std::cout << "Line is " << l << std::endl;
    std::cout << "Sphere is " << table_->sphere(s) << std::endl;
    Sphere_point_3 sa(table_->sphere(s), l);
    if (!sa.is_valid()) return Event_pair();
    Event_point_3 pa(sa);
    Event_point_3 pb(Sphere_point_3(table_->sphere(s), l.opposite()));
    if (pb < pa) std::swap(pa, pb);
    std::cout << "Ret is " << pa << " " << pb << std::endl;
    return Event_pair(pa, pb);
  }
}


CGAL_END_NAMESPACE
