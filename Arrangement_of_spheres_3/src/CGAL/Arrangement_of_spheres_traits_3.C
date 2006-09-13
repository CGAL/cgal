#include <CGAL/Arrangement_of_spheres_traits_3.h>

/* 
   predictes--------------------------------------------------------
*/


bool Arrangement_of_spheres_traits_3::intersects(Key a,
						 Key b) const {
  //Vector_3 d= table_->center(a) - table_->center(b);
  //FT d2= d*d;
  //FT sr2= table_->sphere(a).squared_radius() + table_->sphere(b).squared_radius();
  //return d2 <= sr2;
  /*Vector_3 n=2*(table_->center(a)-table_->center(b));
  FT d= table_->discriminant(b) - table_->discriminant(a);
  FT dda= n*(table_->center(a)-CGAL::ORIGIN);
  FT diff= (d-dda);
  std::cout << n << std::endl << d << " " << dda << " " << diff << std::endl;
  bool ret= CGAL::square(diff) <=  table_->sphere(a).squared_radius(); 
  {*/
  //if (table_->center(a) == table_->center(b)) return false;
  try {
    Point_3 eqpt= table_->equipower_point(a,b);
    std::cout << eqpt  << std::endl;
    CGAL::Bounded_side bs=table_->sphere(a).bounded_side(eqpt);
    CGAL::Bounded_side obs=table_->sphere(b).bounded_side(eqpt);
    return bs != CGAL::ON_UNBOUNDED_SIDE;
  } catch (Sphere_3_table::Equal_centers_exception) {
    std::cout << "Equal centers for " << a << " and " << b << std::endl;
    return false;
  }
  /*}
    return ret;*/
}

bool Arrangement_of_spheres_traits_3::intersects(Key a,
						 Key b,
						 Key c) const {
  CGAL_precondition(intersects(a,b));
  CGAL_precondition(intersects(b,c));
  CGAL_precondition(intersects(c,a));
  // NOTE could be better, I think
  Plane_3 pab = table_->equipower_plane(a,b);
  Plane_3 pac = table_->equipower_plane(a,c);
  CGAL::Object o= di_(pab, pac);
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


bool Arrangement_of_spheres_traits_3::sphere_intersects_sweep(Key sphere,
							      Sphere_point_3 ep) const {
  Event_pair be= sphere_events(sphere);
  if (be.first.compare(ep, sweep_coordinate()) != CGAL::LARGER 
      && be.second.compare(ep, sweep_coordinate()) != CGAL::SMALLER) return true;
  else return false;
}

CGAL::Comparison_result  Arrangement_of_spheres_traits_3::compare_depths(const Sphere_point_3 &a, 
									 const Sphere_point_3 &b) const {
  return a.compare(b, sweep_coordinate());
}

bool Arrangement_of_spheres_traits_3::sphere_intersects_rule(Key s, 
							     Key rule, 
							     Coordinate_index C) const {
  FT d= table_->center(s)[C.index()]-table_->center(rule)[C.index()];
  if (CGAL::square(d) <= table_->sphere(s).squared_radius()) return true;
  else return false;
}


CGAL::Comparison_result 
Arrangement_of_spheres_traits_3::compare_sphere_center_c(Key a,
							 const Sphere_point_3& d,
							 Coordinate_index C) const {
  return CGAL::Comparison_result(-d.compare(table_->center(a), C));
}

CGAL::Comparison_result 
Arrangement_of_spheres_traits_3::compare_sphere_center_c(Key a,
							 FT d,
							 Coordinate_index C) const {
  return CGAL::compare(table_->center(a)[C.index()], d);
}


CGAL::Comparison_result 
Arrangement_of_spheres_traits_3::compare_equipower_point_to_rule(Key a, 
								 Key b,
								 Key c, 
								 Coordinate_index C) const{
  // NOTE redo this with computing it directly
  Point_3 pt= table_->equipower_point(a,b);
  return CGAL::compare(pt[C.index()], table_->center(c)[C.index()]);

}


CGAL::Sign Arrangement_of_spheres_traits_3::sign_of_separating_plane_normal_c(Key a, Key b, 
							       Coordinate_index C) const {
  CGAL_precondition(CGAL::enum_cast<CGAL::Sign>(CGAL::LARGER) == CGAL::POSITIVE);
  int c= compare_sphere_centers_c(b, a, other_plane_coordinate(C));
  if (C==plane_coordinate(0)) c=-c;
  CGAL::Sign ret= CGAL::enum_cast<CGAL::Sign>(c);
  CGAL_assertion(ret== CGAL::sign(table_->separating_plane(a,b).orthogonal_vector()[C.index()]));
  return ret;
  /*if (C==0) {
    sn = CGAL::sign(-table_->center_c(b,1) + table_->center_c(a,1));
  } else {
    sn = CGAL::sign( table_->center_c(b,0) - table_->center_c(a,0));
  }
  return sn;*/
}



CGAL::Sign Arrangement_of_spheres_traits_3::sign_of_equipower_plane_normal_c(Key a, 
						   Key b, Coordinate_index C) const {
  CGAL::Sign sn= CGAL::sign(table_->center(a)[C.index()]-table_->center(b)[C.index()]);
  CGAL_assertion(sn== CGAL::sign(table_->equipower_plane(a,b).orthogonal_vector()[C.index()]));
  return sn;
}



CGAL::Oriented_side Arrangement_of_spheres_traits_3::oriented_side_of_equipower_plane(Key a, Key b,
							    const Sphere_point_3 &s) const {
  CGAL_assertion(s.is_valid());
  Plane_3 p= table_->equipower_plane(a,b);
  CGAL::Object o= di_(p, s.line());
  Point_3 pt;
  Line_3 l;
  if (CGAL::assign(pt, o)){
    CGAL::Comparison_result c= s.compare_on_line(pt);
    FT d= p.orthogonal_vector() * s.line().to_vector();
    if (d > 0) return CGAL::enum_cast<CGAL::Oriented_side>(c);
    else return CGAL::enum_cast<CGAL::Oriented_side>(-c);
  } else if (CGAL::assign( l, o)){
    return CGAL::ON_ORIENTED_BOUNDARY;
  } else {
    return p.oriented_side(s.line().point());
  }
}


CGAL::Oriented_side Arrangement_of_spheres_traits_3::oriented_side_of_center_plane(Key a, Key b,
							 Key sphere_center) const {
  Vector_3 d(table_->center(b)-table_->center(a));
  Line_2 l(Point_2(table_->center(a)[plane_coordinate(0).index()],
		   table_->center(a)[plane_coordinate(1).index()]), 
	   Vector_2(d[plane_coordinate(0).index()], d[plane_coordinate(1).index()]));
  return l.oriented_side(Point_2(table_->center(sphere_center)[plane_coordinate(0).index()],
				 table_->center(sphere_center)[plane_coordinate(1).index()]));
}

CGAL::Comparison_result Arrangement_of_spheres_traits_3::compare_sphere_centers_c(Key a, Key b, Coordinate_index C) const {
  return CGAL::compare(table_->center(a)[C.index()], table_->center(b)[C.index()]);
  }





CGAL::Bounded_side 
Arrangement_of_spheres_traits_3::bounded_side_of_sphere_on_equipower_plane_rule_line(Key sphere0,
										     Key sphere1,
										     Key rule,
										     Coordinate_index C,
										     const Sphere_point_3& sp) const {
  
  Plane_3 eqp= table_->equipower_plane(sphere0, sphere1);
  FT v[3], p[3];
  v[C.index()]=1;
  v[other_plane_coordinate(C).index()]=0;
  v[sweep_coordinate().index()] = 0;
  p[plane_coordinate(0).index()]=table_->center(rule)[plane_coordinate(0).index()];
  p[plane_coordinate(1).index()]=table_->center(rule)[plane_coordinate(1).index()];
  p[sweep_coordinate().index()]= 0;
  Plane_3 pp(Point_3(p[0],p[1],p[2]),
	     Vector_3(v[0], v[1], v[2]));
  Line_3 l;
  CGAL::Object o= di_(eqp, pp);
  CGAL_assertion(CGAL::assign(l,o));
  CGAL::assign(l,o);
  std::cout << "The intersection line is " << l << std::endl;
  //CGAL::Comparison_result dir= CGAL::compare(ep.simple_coordinate(1-C), spheres_[h->curve().key()].table_->center(1-C));*/
  Sphere_point_3 p0(table_->sphere(sphere0), l);


  if (!p0.is_valid()){
    return CGAL::ON_UNBOUNDED_SIDE;
  } else {
    Sphere_point_3 p1(table_->sphere(sphere0), l.opposite());
    
    if (compare_depths(p0, p1)== CGAL::LARGER) {
      std::swap(p0, p1);
    }
    //CGAL_assertion(p0 < p1);
    if (compare_depths(sp,p0)==CGAL::EQUAL 
	|| compare_depths(sp, p1) == CGAL::EQUAL) {
      return CGAL::ON_BOUNDARY;
    }
    if (compare_depths(sp, p0) == CGAL::LARGER 
	&& compare_depths(sp, p1) == CGAL::SMALLER) return CGAL::ON_BOUNDED_SIDE;
    else return CGAL::ON_UNBOUNDED_SIDE;
  }
}


  // compare along C where the intersection point of the two rules
  // (defined by Coord C of ep) is relative to the sphere of the arc
CGAL::Bounded_side 
Arrangement_of_spheres_traits_3::bounded_side_of_sphere(Key s,
							Key x,
							Key y,
							const Sphere_point_3 &z) const{
  FT p[3];
  p[plane_coordinate(0).index()]= table_->center(x)[plane_coordinate(0).index()];
  p[plane_coordinate(1).index()]= table_->center(y)[plane_coordinate(1).index()];
  p[sweep_coordinate().index()]=0;
  FT v[3];
  v[plane_coordinate(0).index()]=0;
  v[plane_coordinate(1).index()]=0;
  v[sweep_coordinate().index()]=1;
  /*std::cout << "Line through " << CGAL::to_double(p[0]) 
    << " " << CGAL::to_double(p[1]) << std::endl;*/
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





/*
  Quadratic constructions--------------------------------------------------
*/

Arrangement_of_spheres_traits_3::Quadratic_NT 
Arrangement_of_spheres_traits_3::intersection_c(Key s, 
						Line_3 l, 
						Coordinate_index C) const {
  //sort on correct coordinate
  bool first= CGAL::sign(l.to_vector()[C.index()]) == CGAL::POSITIVE;

  Vector_3 lp= l.point()-CGAL::ORIGIN;
  Vector_3 vc= table_->center(s)-CGAL::ORIGIN; 
  Vector_3 lv= l.to_vector(); 
  FT a=lv*lv;
  FT b=2*lv*(lp-vc); //-2*lv* vc + 2*lv*lp;
  FT c=lp*lp + vc*vc-table_->sphere(s).squared_radius()-2*lp*vc;

  FT disc= b*b-4*a*c;
  CGAL_assertion(disc >= 0);
  CGAL_precondition(a!=0);

  bool rt= first && CGAL::sign(lv[C.index()]) == CGAL::POSITIVE 
    || !first && CGAL::sign(lv[C.index()]) == CGAL::NEGATIVE;
  Quadratic_NT t=CGAL::make_root_of_2(a,b,c,rt);

   
  if (first) {
    CGAL_assertion(t*lv[C.index()] <= lv[C.index()] *CGAL::make_root_of_2(a,b,c, !rt));
  } else {
    CGAL_assertion(t*lv[C.index()] >= lv[C.index()] *CGAL::make_root_of_2(a,b,c, !rt)); 
  }
  Quadratic_NT p0= lp[C.index()] + lv[C.index()]*t;
  //std::cout << p0 << " is " << CGAL::to_double(p0) << std::endl;
  return p0;
}



Arrangement_of_spheres_traits_3::Event_pair 
Arrangement_of_spheres_traits_3::sphere_events(Key ind) const {
  FT v[3];
  v[plane_coordinate(0).index()]=0;
  v[plane_coordinate(1).index()]=0;
  v[sweep_coordinate().index()]=1;
  return Event_pair(Sphere_point_3(table_->sphere(ind), 
				   Line_3(table_->center(ind), Vector_3(v[0], v[1], v[2]))),
		    Sphere_point_3(table_->sphere(ind), 
				   Line_3(table_->center(ind), Vector_3(-v[0], -v[1], -v[2]))));
}


Arrangement_of_spheres_traits_3::Event_pair 
Arrangement_of_spheres_traits_3::intersection_2_events(Key a, Key b) const {
  std::cout << "Intersection events for " << a << " and " << b << std::endl;
  if (!intersects(a,b)) {
    return Event_pair();
  }
  Plane_3 eqp= table_->equipower_plane(a,b);
  Plane_3 sp= table_->separating_plane(a,b);
  //CGAL_assertion(eqp.has_on(eqpoint));
  CGAL::Object o= di_(eqp, sp);
  std::cout << "eqp= " << eqp << std::endl << "sp= " << sp << std::endl;
    //<< "eqpoint= " << eqpoint << std::endl;
  Line_3 l;
  if (!CGAL::assign(l,o)){
    CGAL_assertion(0);
  }
  Line_3 lf, lb;
  
  if (l.to_vector()[sweep_coordinate().index()] > 0){
    lf=l;
    lb=l.opposite();
  } else {
    lf=l.opposite();
    lb=l;
  }
  if (l.to_vector()[sweep_coordinate().index()]==0) {
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
  CGAL_assertion(compare_depths(ret.first, ret.second) != CGAL::LARGER);
  return ret;
}

Arrangement_of_spheres_traits_3::Event_pair 
Arrangement_of_spheres_traits_3::intersection_3_events(Key a, Key b, Key c) const {
  CGAL_precondition(a != b && a != c && c!= a);
  if (!intersects(a,b, c)) {
    return Event_pair();
  }
  Plane_3 eqpab= table_->equipower_plane(a,b);
  Plane_3 eqpac= table_->equipower_plane(a,c);

  CGAL::Object o= di_(eqpab, eqpac);
  Line_3 l;
  if (!CGAL::assign(l,o)){
    // they intersect on the same circle. Ick.
    CGAL_assertion(0);
  }
  Line_3 lf, lb;
  CGAL_assertion(l.to_vector()[sweep_coordinate().index()] != 0);
  if (l.to_vector()[sweep_coordinate().index()] > 0){
    lf=l;
    lb=l.opposite();
  } else {
    lf=l.opposite();
    lb=l;
  }
  if (l.to_vector()[sweep_coordinate().index()] ==0) {
    std::cerr << "Degeneracy, picking one for start of threesome " 
	      << a << " " << b << " " << c << std::endl;
  }
  
  Event_pair ret(Sphere_point_3(table_->sphere(a), lf),
		 Sphere_point_3(table_->sphere(a), lb));
  CGAL_assertion(compare_depths(ret.first, ret.second) != CGAL::LARGER);
  return ret;
}

Arrangement_of_spheres_traits_3::Event_pair 
Arrangement_of_spheres_traits_3::sphere_intersect_extremum_events(Key a, 
								  Coordinate_index C,
								  Key b) const {
  CGAL_precondition(intersects(a,b));
  Plane_3 eqpab= table_->equipower_plane(a,b);
  Plane_3 rp= rule_plane(a,C);
  CGAL::Object o= di_(eqpab, rp);
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

Arrangement_of_spheres_traits_3::Event_pair 
Arrangement_of_spheres_traits_3::sphere_intersect_rule_rule_events(Key s, 
								   Key x,
								   Key y) const {
  CGAL_precondition(s != Key());
  CGAL_precondition(x != Key());
  CGAL_precondition(y != Key());

  std::cout << "SRR event for " << s << " " << x << " " << y << std::endl;

  Plane_3 x_plane= rule_plane(x,plane_coordinate(0));
  Plane_3 y_plane= rule_plane(y,plane_coordinate(1));
  CGAL::Object o= di_(x_plane, y_plane);
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

Arrangement_of_spheres_traits_3::Event_pair 
Arrangement_of_spheres_traits_3::sphere_intersect_rule_events(Key a, Key r, Coordinate_index C) const {
  CGAL_precondition(sphere_intersects_rule(a, r, C));

  FT pa[3];
  pa[sweep_coordinate().index()]= table_->center(a)[sweep_coordinate().index()];
  pa[C.index()]=  table_->center(r)[C.index()];
  pa[other_plane_coordinate(C).index()]
    =  table_->center(a)[other_plane_coordinate(C).index()];

  FT va[3]={0,0,0};
  va[sweep_coordinate().index()]= 1;
  Line_3 l(Point_3(pa[0], pa[1], pa[2]),
	   Vector_3(va[0], va[1], va[2]));
  Sphere_point_3 sp(table_->sphere(a), l);
  Sphere_point_3 spn(table_->sphere(a), l.opposite());
  return Event_pair(sp, spn);
}


Arrangement_of_spheres_traits_3::Event_pair 
Arrangement_of_spheres_traits_3::circle_cross_rule_events(Key a, Key b,
							  Key rs, Coordinate_index C) const {
  Plane_3 eqp= table_->equipower_plane(a,b);
  /*FT na[3]={0,0,0};
    na[C.index()]=1;*/
  Plane_3 rp = rule_plane(rs, C); //(table_->center(rs), Vector_3(na[0], na[1], na[2]));
  CGAL::Object o=di_(eqp, rp);
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


Arrangement_of_spheres_traits_3::Plane_3
Arrangement_of_spheres_traits_3::rule_plane(Key r, Coordinate_index C) const {
  FT n[3]={0,0,0};
  n[C.index()]= 1;
  return Plane_3(table_->center(r), Vector_3(n[0], n[1], n[2]));
}
