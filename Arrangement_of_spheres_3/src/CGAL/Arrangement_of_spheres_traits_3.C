#include <CGAL/Arrangement_of_spheres_traits_3.h>



void Arrangement_of_spheres_traits_3::initialize_1() {
  di_ = geometric_traits_object().intersect_3_object();
  spheres_.push_back(Sphere_3());
  spheres_.push_back(Sphere_3());
  spheres_.push_back(Sphere_3());
}


void Arrangement_of_spheres_traits_3::initialize_2() {
   CGAL::Bbox_3 b(std::numeric_limits<double>::max(),
		  std::numeric_limits<double>::max(),
		  std::numeric_limits<double>::max(),
		  -std::numeric_limits<double>::max(),
		  -std::numeric_limits<double>::max(),
		  -std::numeric_limits<double>::max());
   for (unsigned int i=3; i< spheres_.size(); ++i){
     b= b+ spheres_[i].bbox();
   }
   FT inf=2*std::max(b.xmax(),
		     std::max(std::abs(b.xmin()),
			      std::max(b.ymax(),
				       std::max(std::abs(b.ymin()),
						std::max(b.zmax(),
							 std::abs(b.zmin()))))));
   spheres_[1]=Sphere_3(Point_3(-inf, -inf, -inf), 0);
   spheres_[2]=Sphere_3(Point_3(inf, inf, inf), 0);
}

/* 
   predictes--------------------------------------------------------
*/


bool Arrangement_of_spheres_traits_3::intersects(Key a,
						 Key b) const {
  Vector_3 d= center(a) - center(b);
  FT d2= d*d;
  FT sr2= sphere(a).squared_radius() + sphere(b).squared_radius();
  return d2 <= sr2;
}

bool Arrangement_of_spheres_traits_3::intersects(Key a,
						 Key b,
						 Key c) const {
  CGAL_precondition(intersects(a,b));
  CGAL_precondition(intersects(b,c));
  CGAL_precondition(intersects(c,a));
  // NOTE could be better, I think
  Plane_3 pab = equipower_plane(a,b);
  Plane_3 pac = equipower_plane(a,c);
  CGAL::Object o= di_(pab, pac);
  Line_3 l;
  if (CGAL::assign(l,o)){
    Sphere_point_3 sp(sphere(a), l);
    return sp.is_valid();
  } else {
    CGAL_assertion(pab== pac);
    return true;
  }
  
}


CGAL::Comparison_result  Arrangement_of_spheres_traits_3::compare_depths(const Sphere_point_3 &a, 
									 const Sphere_point_3 &b) const {
  return a.compare(b, Coordinate_index(2));
}

bool Arrangement_of_spheres_traits_3::sphere_intersects_rule(Key s, Key rule, Coordinate_index C) const {
  FT d= center_c(s, C)-center_c(rule,C);
  if (CGAL::square(d) <= sphere(s).squared_radius()) return true;
  else return false;
}



CGAL::Comparison_result 
Arrangement_of_spheres_traits_3::compare_equipower_point_to_rule(Key a, 
								 Key b,
								 Key c, 
								 Coordinate_index C) const{
  // NOTE redo this with computing it directly
  Point_3 pt= equipower_point(a,b);
  return CGAL::compare(pt[C], center_c(c, C));

}


CGAL::Sign Arrangement_of_spheres_traits_3::sign_of_separating_plane_normal_c(Key a, Key b, 
							       Coordinate_index C) const {
  CGAL_precondition(CGAL::enum_cast<CGAL::Sign>(CGAL::LARGER) == CGAL::POSITIVE);
  int c= compare_sphere_centers_c(b, a,1-C);
  if (C==0) c=-c;
  CGAL::Sign ret= CGAL::enum_cast<CGAL::Sign>(c);
  CGAL_assertion(ret== CGAL::sign(separating_plane(a,b).orthogonal_vector()[C]));
  return ret;
  /*if (C==0) {
    sn = CGAL::sign(-center_c(b,1) + center_c(a,1));
  } else {
    sn = CGAL::sign( center_c(b,0) - center_c(a,0));
  }
  return sn;*/
}



CGAL::Sign Arrangement_of_spheres_traits_3::sign_of_equipower_plane_normal_c(Key a, 
						   Key b, Coordinate_index C) const {
  CGAL::Sign sn= CGAL::sign(center_c(a,C)-center_c(b,C));
  CGAL_assertion(sn== CGAL::sign(equipower_plane(a,b).orthogonal_vector()[C]));
  return sn;
}



CGAL::Oriented_side Arrangement_of_spheres_traits_3::oriented_side_of_equipower_plane(Key a, Key b,
							    const Sphere_point_3 &s) const {
  CGAL_assertion(s.is_valid());
  Plane_3 p= equipower_plane(a,b);
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
  Vector_3 d(center(b)-center(a));
  Line_2 l(Point_2(center_c(a,Coordinate_index(0)),
			 center_c(a,Coordinate_index(1))), 
	      Vector_2(d.x(), d.y()));
  return l.oriented_side(Point_2(center_c(sphere_center, Coordinate_index(0)),
				    center_c(sphere_center, Coordinate_index(1))));
}

CGAL::Comparison_result Arrangement_of_spheres_traits_3::compare_sphere_centers_c(Key a, Key b, Coordinate_index C) const {
  return CGAL::compare(center_c(a,C), center_c(b,C));
  }





CGAL::Bounded_side Arrangement_of_spheres_traits_3::bounded_side_of_sphere_on_equipower_plane_rule_line(Key sphere0,
													Key sphere1,
													Key rule,
													Coordinate_index C,
													const Sphere_point_3& sp) const {
  
  Plane_3 eqp= equipower_plane(sphere0, sphere1);
  FT v[2];
  v[C]=1;
  v[1-C]=0;
  Plane_3 pp(Point_3(center_c(rule,Coordinate_index(0)),
		     center_c(rule,Coordinate_index(1)), 0),
	     Vector_3(v[0], v[1], 0));
  Line_3 l;
  CGAL::Object o= di_(eqp, pp);
  CGAL_assertion(CGAL::assign(l,o));
  CGAL::assign(l,o);
  std::cout << "The intersection line is " << l << std::endl;
  //CGAL::Comparison_result dir= CGAL::compare(ep.simple_coordinate(1-C), spheres_[h->curve().key()].center(1-C));*/
  Sphere_point_3 p0(sphere(sphere0), l);


  if (!p0.is_valid()){
    return CGAL::ON_UNBOUNDED_SIDE;
  } else {
    Sphere_point_3 p1(sphere(sphere0), l.opposite());
    
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
  FT p[2];
  p[0]= center_c(x, Coordinate_index(0));
  p[1]= center_c(y, Coordinate_index(1));
    
  /*std::cout << "Line through " << CGAL::to_double(p[0]) 
    << " " << CGAL::to_double(p[1]) << std::endl;*/
  Line_3 l(Point_3(p[0], p[1], 0), Vector_3(0,0,1));
    
  Sphere_point_3 p0(sphere(s), l);

  if (!p0.is_valid()) {
    return CGAL::ON_UNBOUNDED_SIDE;
  } else {
    Sphere_point_3 p1(sphere(s), 
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
   Constructions--------------------------------------------------------
*/


Arrangement_of_spheres_traits_3::Quadratic_NT 
Arrangement_of_spheres_traits_3::intersection_c(Key s, 
						Line_3 l, 
						Coordinate_index C) const {
  //sort on correct coordinate
  bool first= CGAL::sign(l.to_vector()[C]) == CGAL::POSITIVE;

  Vector_3 lp= l.point()-CGAL::ORIGIN;
  Vector_3 vc= center(s)-CGAL::ORIGIN; 
  Vector_3 lv= l.to_vector(); 
  FT a=lv*lv;
  FT b=2*lv*(lp-vc); //-2*lv* vc + 2*lv*lp;
  FT c=lp*lp + vc*vc-sphere(s).squared_radius()-2*lp*vc;

  FT disc= b*b-4*a*c;
  CGAL_assertion(disc >= 0);
  CGAL_precondition(a!=0);

  bool rt= first && CGAL::sign(lv[C]) == CGAL::POSITIVE 
    || !first && CGAL::sign(lv[C]) == CGAL::NEGATIVE;
  Quadratic_NT t=CGAL::make_root_of_2(a,b,c,rt);

   
  if (first) {
    CGAL_assertion(t*lv[C] <= lv[C] *CGAL::make_root_of_2(a,b,c, !rt));
  } else {
    CGAL_assertion(t*lv[C] >= lv[C] *CGAL::make_root_of_2(a,b,c, !rt)); 
  }
  Quadratic_NT p0= lp[C] + lv[C]*t;
  //std::cout << p0 << " is " << CGAL::to_double(p0) << std::endl;
  return p0;
}


Arrangement_of_spheres_traits_3::Point_3 
Arrangement_of_spheres_traits_3::center(Key a) const {
  return sphere(a).center();
}


Arrangement_of_spheres_traits_3::Plane_3 
Arrangement_of_spheres_traits_3::separating_plane(Key a, Key b) const {
  Vector_3 d(center(b)-center(a));
  Vector_3 n(-d[1], d[0], 0);
  return Plane_3(center(b), n);
}



Arrangement_of_spheres_traits_3::Plane_3 
Arrangement_of_spheres_traits_3::equipower_plane(Key a, Key b) const {
  CGAL_precondition(a != b);
  Vector_3 n=2*(center(a)-center(b));
  FT d= disc(b) - disc(a);
  return Plane_3(n[0], n[1], n[2], d);
}

Arrangement_of_spheres_traits_3::Point_3 
Arrangement_of_spheres_traits_3::equipower_point(Key a, Key b) const {
  Plane_3 eqp= equipower_plane(a,b);
  Line_3 l(center(a), (center(a)-center(b)));
  CGAL::Object o= di_(eqp, l);
  Point_3 pt;
  CGAL_assertion(CGAL::assign(pt, o));
  CGAL::assign(pt, o);
  return pt;
}

Arrangement_of_spheres_traits_3::FT 
Arrangement_of_spheres_traits_3::disc(Key i) const {
  Vector_3 v= center(i)-CGAL::ORIGIN;
  return v*v - sphere(i).squared_radius();
}
  

Arrangement_of_spheres_traits_3::FT 
Arrangement_of_spheres_traits_3::center_c(Key sphere,
						Coordinate_index C) const {
  //CGAL_precondition(c.is_rule());
  return center(sphere)[C];
}




Arrangement_of_spheres_traits_3::Sphere_3 Arrangement_of_spheres_traits_3::sphere(Key ind) const {
  //  CGAL_precondition(ind+3 >=0);
  //CGAL_precondition(ind+3 < spheres_.size());
  //CGAL_precondition(ind != TEMP_ID || has_temp_);
  return spheres_[ind.internal_index()];
  /*if (a==BL_ID) {
    return bl_sphere_;
  } else if (a==TR_ID) {
    return tr_sphere_;
  } else if (a== TEMP_ID) {
    CGAL_assertion(has_temp_);
    return temp_sphere_;
  } else {
    CGAL_precondition(a < spheres_.size());
    CGAL_precondition(a >= 0);
    return spheres_[a];
    }*/
}





/*
  Quadratic constructions--------------------------------------------------
*/





Arrangement_of_spheres_traits_3::Event_pair 
Arrangement_of_spheres_traits_3::sphere_events(Key ind) const {
  return Event_pair(Sphere_point_3(sphere(ind), 
				   Line_3(center(ind), Vector_3(0,0,1))),
		    Sphere_point_3(sphere(ind), 
				   Line_3(center(ind), Vector_3(0,0,-1))));
}


Arrangement_of_spheres_traits_3::Event_pair 
Arrangement_of_spheres_traits_3::intersection_events(Key a, Key b) const {
  CGAL_precondition(intersects(a,b));
  Plane_3 eqp= equipower_plane(a,b);
  Plane_3 fp(Point_3(equipower_point(a,b)), Vector_3(1,0,0));
  CGAL::Object o= di_(eqp, fp);
  Line_3 l;
  if (!CGAL::assign(l,o)){
    Plane_3 fp2(Point_3(equipower_point(a,b)), Vector_3(0,1,0));
    CGAL::Object o= di_(eqp, fp);
    if (!CGAL::assign(l,o)){
      CGAL_assertion(0);
    }
  }
  Line_3 lf, lb;
  
  if (l.to_vector().z() > 0){
    lf=l;
    lb=l.opposite();
  } else {
    lf=l.opposite();
    lb=l;
  }
  if (l.to_vector().z()==0) {
    std::cerr << "Degeneracy, picking one for start of circle " 
	      << a << " " << b << std::endl;
  }
  
  Event_pair ret(Sphere_point_3(sphere(a), lf),
		 Sphere_point_3(sphere(a), lb));
  CGAL_assertion(compare_depths(ret.first, ret.second) != CGAL::LARGER);
  return ret;
}

Arrangement_of_spheres_traits_3::Event_pair 
Arrangement_of_spheres_traits_3::intersection_events(Key a, Key b, Key c) const {
  CGAL_precondition(intersects(a,b, c));
  Plane_3 eqpab= equipower_plane(a,b);
  Plane_3 eqpac= equipower_plane(a,c);

  CGAL::Object o= di_(eqpab, eqpac);
  Line_3 l;
  if (!CGAL::assign(l,o)){
    // they intersect on the same circle. Ick.
    CGAL_assertion(0);
  }
  Line_3 lf, lb;
  CGAL_assertion(l.to_vector().z() != 0);
  if (l.to_vector().z() > 0){
    lf=l;
    lb=l.opposite();
  } else {
    lf=l.opposite();
    lb=l;
  }
  if (l.to_vector().z()==0) {
    std::cerr << "Degeneracy, picking one for start of threesome " 
	      << a << " " << b << " " << c << std::endl;
  }
  
  Event_pair ret(Sphere_point_3(sphere(a), lf),
		 Sphere_point_3(sphere(a), lb));
  CGAL_assertion(compare_depths(ret.first, ret.second) != CGAL::LARGER);
  return ret;
}
