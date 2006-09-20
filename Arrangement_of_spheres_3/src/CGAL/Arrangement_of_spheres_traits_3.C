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
    //CGAL::Bounded_side obs=table_->sphere(b).bounded_side(eqpt);
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
							      const Sphere_point_3& ep) const {
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


bool Arrangement_of_spheres_traits_3::sphere_intersects_rule(Key s, 
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
								 const Sphere_point_3 &sp, 
								 Coordinate_index C) const{
  // NOTE redo this with computing it directly
  Point_3 pt= table_->equipower_point(a,b);
  return CGAL::Comparison_result(-sp.compare(pt, C));
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

CGAL::Oriented_side 
Arrangement_of_spheres_traits_3::oriented_side(const Plane_3 &p,
					       const Sphere_point_3 &s) const{
  CGAL::Object o= di_(p, s.line());
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

CGAL::Oriented_side Arrangement_of_spheres_traits_3::oriented_side_of_equipower_plane(Key a, Key b,
							    const Sphere_point_3 &s) const {
  CGAL_assertion(s.is_valid());
  Plane_3 p= table_->equipower_plane(a,b);
  return oriented_side(p, s);
}


CGAL::Oriented_side 
Arrangement_of_spheres_traits_3::oriented_side_of_separating_plane(Key a, Key b,
							       const Sphere_point_3 &ep) const {
  Plane_3 sp= table_->separating_plane(a,b);
  return oriented_side(sp, ep);
}

CGAL::Comparison_result Arrangement_of_spheres_traits_3::compare_sphere_centers_c(Key a, Key b, Coordinate_index C) const {
  return CGAL::compare(table_->center(a)[C.index()], table_->center(b)[C.index()]);
  }



CGAL::Comparison_result 
Arrangement_of_spheres_traits_3::compare_sphere_sphere_at_sweep(const Sphere_point_3 &t,
								Key sphere0,
								Key sphere1,
								const Sphere_point_3 &pt,
								Coordinate_index C) const {
   {
     std::cout << "Compare SS at sweep " << t << " spheres " 
	       <<  sphere0 << " " << sphere1 
	       << " point " << pt << " line " << pt.line() << " coord " << C 
	       << std::endl;

     bool i0= sphere_intersects_rule(sphere0, pt, C);
     bool i1= sphere_intersects_rule(sphere1, pt, C);
     
     
     //std::cout << "Comparisons are " << hi << " " << hp << " and " << ohi << " " << ohp << std::endl;
     if (!i0) {
       CGAL::Comparison_result sc0= compare_sphere_center_c(sphere0, pt, C);
       std::cout << "Returning " << sc0 << " since sphere 0 misses" << std::endl;
       return sc0;
     } else if (!i1) {
       CGAL::Comparison_result sc1= compare_sphere_center_c(sphere1, pt, C);
       std::cout << "Returning " << sc1 << " since sphere 0 misses" << std::endl;
       return sc1;
     }
   }


   Plane_3 eqp= table_->equipower_plane(sphere0, sphere1);
   std::cout << "EQP is " << eqp << std::endl;
   
   FT v[3]={0,0,0};
   v[other_plane_coordinate(C).index()]=1;
   
   Plane_3 rule_plane(pt.line().point(),
		      pt.line().point()+ Vector_3(v[0], v[1], v[2]),
		      pt.line().point()+ pt.line().to_vector());
   if (rule_plane.orthogonal_vector()[C.index()] < 0) {
     rule_plane=rule_plane.opposite();
   }
   std::cout << "rule_plane is " << rule_plane << std::endl;
   CGAL_assertion(rule_plane != Plane_3());
   Line_3 l;
   CGAL::Object o= di_(eqp, rule_plane);
   if (!CGAL::assign(l,o)) {
     CGAL::Comparison_result cr= compare_equipower_point_to_rule(sphere0, sphere1,
								 pt, C);
     std::cout << "Returning " << cr << " since sphere the planes miss" << std::endl;
     return cr;
   }
 
   std::cout << "The intersection line is " << l << std::endl;
  

   
   Sphere_point_3 switch_point_0(table_->sphere(sphere0), l);
   
   if (!switch_point_0.is_valid()){
     
     CGAL::Comparison_result cr=  compare_equipower_point_to_rule(sphere0,sphere1,
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
   CGAL::Comparison_result spcr0= switch_point_0.compare(t, sweep_coordinate());
   
   Sphere_point_3 switch_point_1(table_->sphere(sphere0), l.opposite());
   std::cout << "switch points are " << switch_point_0 
	     << " and " << switch_point_1 << std::endl;
   CGAL::Comparison_result spcr1= switch_point_1.compare(t, sweep_coordinate());
   std::cout << "spcr is " << spcr0 << spcr1 << std::endl;
   
   
   CGAL::Oriented_side os_0= oriented_side(sepp, switch_point_0);
   CGAL::Oriented_side os_1= oriented_side(sepp, switch_point_1);
   std::cout << "os are " << os_0 << os_1 << std::endl;
   
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




CGAL::Bounded_side 
Arrangement_of_spheres_traits_3
::bounded_side_of_sphere_projected(const Sphere_point_3 &t,
				   Key sphere,
				   Key oplane,
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
  CGAL::Object o= di_(rp, p);
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

  // compare along C where the intersection point of the two rules
  // (defined by Coord C of ep) is relative to the sphere of the arc
CGAL::Bounded_side 
Arrangement_of_spheres_traits_3::bounded_side_of_sphere(Key s,
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



bool
Arrangement_of_spheres_traits_3::is_over_circle_c(Key s,
						  const Sphere_point_3 &z,
						  Coordinate_index C) const {
  Line_3 l= z.line();
  FT p[3], v[3];
  for (unsigned int i=0; i< 3; ++i) {
    p[i]= l.point()[i];
    v[i]= l.to_vector()[i];
  }
  Coordinate_index oc= other_plane_coordinate(C);
  p[oc.index()]= table_->center(s)[oc.index()];
  v[oc.index()]= 0;
  Line_3 nl=Line_3(Point_3(p[0], p[1], p[2]),Vector_3(v[0], v[1], v[2]));
  Sphere_point_3 sp(table_->sphere(s), nl);
  if (!sp.is_valid()) return false;
  else {
    Sphere_point_3 spo(table_->sphere(s), nl.opposite());
    CGAL::Comparison_result c= sp.compare(z, sweep_coordinate());
    CGAL::Comparison_result co= spo.compare(z, sweep_coordinate());
    if (co == c) return true;
    else return false;
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
