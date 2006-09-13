#include <CGAL/Arrangement_of_spheres_3/Slice.h>

#define DPRINT(x) x





/* 
   2D -------------------------------------------------------------
*/

bool Slice::intersects_rz(T::Key a, NT z) const {
  T::FT p[3];
  p[plane_coordinate(0).index()]=0;
  p[plane_coordinate(1).index()]=1;
  p[sweep_coordinate().index()]=z;
  return t_.sphere_intersects_sweep(a, T::Sphere_point_3(T::Point_3(p[0], p[1], p[2]))); 
}


Slice::T::Circle_2 Slice::circle_rz(T::Key a, NT z) const {
  T::Sphere_3 s= t_.sphere(a);
  NT r2= t_.sphere(a).squared_radius()
    -  CGAL::square(t_.sphere(a).center()[sweep_coordinate().index()]-z);
  if (r2 >=0) {
    T::Circle_2  c(T::Point_2(t_.sphere(a).center()[plane_coordinate(0).index()], 
			      t_.sphere(a).center()[plane_coordinate(1).index()]), r2);
    return c;
  } else {
    return T::Circle_2(T::Point_2(t_.sphere(a).center()[plane_coordinate(0).index()], 
				  t_.sphere(a).center()[plane_coordinate(1).index()]), 0);
  }
}


Slice::T::Point_2 Slice::center_point_rz(T::Key a, T::Key b, NT z) const {
  //T::Intersect_with_sweep is= tr_.intersect_with_sweep_object(z);
  T::Circle_2 ca= circle_rz(a,z);
  T::Circle_2 cb= circle_rz(b,z);
  T::Geometric_traits::Compute_squared_length_2 csl= t_.geometric_traits_object().compute_squared_length_2_object();
  T::FT c02= csl(ca.center()-CGAL::ORIGIN);
  T::FT c12= csl(cb.center()-CGAL::ORIGIN);
  T::FT c0c1= (cb.center()-CGAL::ORIGIN)*(ca.center()-CGAL::ORIGIN);
  NT sqr=c02 + c12 -2*c0c1;
  NT t=(cb.squared_radius() - ca.squared_radius() +sqr)
    /(2*sqr);
  T::Point_2 p= CGAL::ORIGIN+ (t*(ca.center()-CGAL::ORIGIN) + (1-t)*(cb.center()-CGAL::ORIGIN));
    
  return p;
}



Slice::T::Sphere_point_3 Slice::sphere_point_rz(Sds::Point pt, NT z) const {
  if (pt.is_rule_rule()){
    T::Point_2 p= compute_rule_rule_intersection(pt.rule_key(plane_coordinate(0)),
						 pt.rule_key(plane_coordinate(1)));
    T::FT pt[3]={0,0,0};
    pt[plane_coordinate(0).index()]=p.x();
    pt[plane_coordinate(1).index()]=p.y();
    T::Point_3 p3(pt[0], pt[1], pt[2]);
    T::FT v[3];
    v[plane_coordinate(0).index()]=0;
    v[plane_coordinate(1).index()]=0;
    v[sweep_coordinate().index()]=1;
    T::Sphere_point_3 ret(p3, T::Line_3(p3,
					T::Vector_3(v[0], v[1], v[2])));
    if (!ret.is_valid()) {
      throw p3;
    } else return ret;
  } else if (pt.is_sphere_rule()) {
    T::Sphere_3 s(t_.sphere(pt.sphere_key()));
    T::Line_3 l;
    if (pt.is_smaller()) {
      l= positive_line_rz(pt.rule_key(), pt.rule_coordinate(), z);
    } else {
      l = negative_line_rz(pt.rule_key(), pt.rule_coordinate(), z);
    }
    T::Sphere_point_3 sp (s,l);
    if (!sp.is_valid()) {
      std::cerr << "Constructing point for " << pt << " at " << z << std::endl;
      std::cerr << "Line " << l << " does not intersect sphere " << s << std::endl;
      throw l.projection(s.center());
    } else return sp;
  } else if (pt.sphere_key(0) == pt.sphere_key(1)) {
    CGAL_assertion(0);
    /*int ipt=static_cast<int>(pt.sphere(0).part() & pt.sphere(1).part()) 
      & (~Sds::Curve::ARC_BIT);
    Sds::Curve::Part cpt= static_cast<Sds::Curve::Part>(ipt);
    Sds::Curve rule(pt.sphere(0).key(), cpt);
    CGAL_assertion(rule.is_rule());
    return T::Sphere_point_3(t_.sphere(pt.sphere(0).key()), 
    in_line_rz(rule, z));*/
    return T::Sphere_point_3();
  } else {
    std::cout << "Computing point for " << pt << std::endl;
    CGAL_precondition(pt.sphere_key(0) != pt.sphere_key(1));
    T::Point_2 cp= center_point_rz(pt.sphere_key(0), pt.sphere_key(1), z);
    T::Vector_3 v= t_.sphere(pt.sphere_key(0)).center() 
      - t_.sphere(pt.sphere_key(1)).center();
    std::cout << "cp = " << cp << std::endl;
    std::cout << "v = " << v << std::endl;
    T::FT p[3],vv[3];
    p[plane_coordinate(0).index()]= cp[0];
    p[plane_coordinate(1).index()]= cp[1];
    p[sweep_coordinate().index()]= z;
    vv[plane_coordinate(0).index()]= -v[plane_coordinate(1).index()];
    vv[plane_coordinate(1).index()]= v[plane_coordinate(0).index()];
    vv[sweep_coordinate().index()]=0;
    std::cout << p[0] << " " << p[1] << " " << p[2] << std::endl;
    std::cout << vv[0] << " " << vv[1] << " " << vv[2] << std::endl;
    T::Line_3 l(T::Point_3(p[0], p[1], p[2]), T::Vector_3(vv[0], vv[1], vv[2]));
    
    T::Sphere_point_3 sli(t_.sphere(pt.sphere_key(1)), l);
    if (!sli.is_valid()) {
      std::cerr << "Constructing point for " << pt << " at " << z << std::endl;
      std::cerr << "Line " << l << " does not intersect sphere " 
		<< t_.sphere(pt.sphere_key(1)) << std::endl;
      throw l.projection(t_.sphere(pt.sphere_key(1)).center());
    }
    //std::cout << sli << std::endl;
    //std::cout << sli.approximate_coordinate(0) << ", " <<  sli.approximate_coordinate(1) << std::endl;
    T::Sphere_point_3 osli(t_.sphere(pt.sphere_key(1)), l);
    //std::cout << osli.approximate_coordinate(0) << ", " <<  osli.approximate_coordinate(1) << std::endl;
    
    CGAL_exactness_assertion(sli.compare(osli,Coordinate_index::X()) == CGAL::EQUAL);
    CGAL_exactness_assertion(sli.compare(osli,Coordinate_index::Y()) == CGAL::EQUAL);
    CGAL_exactness_assertion(sli.compare(osli,Coordinate_index::Z()) == CGAL::EQUAL);
    CGAL_assertion(sli.is_valid());
    CGAL_exactness_assertion(CGAL::abs(sli.exact_coordinate(plane_coordinate(0))) 
			     < t_.max_coordinate());
    CGAL_exactness_assertion(CGAL::abs(sli.exact_coordinate(plane_coordinate(1)))
			     < t_.max_coordinate());
    return sli;
  }
}

/* 
   constructions----------------------------------------------------
*/






Slice::T::Point_2 Slice::compute_rule_rule_intersection(T::Key ra,
							T::Key rb) const {
  //CGAL_precondition(ra.is_rule() && rb.is_rule());
  //CGAL_precondition(ra.is_vertical());
  //CGAL_precondition(!rb.is_vertical());
  /*NT x,y;
    if (ra.is_finite()){
    x= spheres_[ra.key()].center().x();
    } else {
    if (!ra.is_negative()) x= inf_;
    else x=-inf_;
    }
    if (rb.is_finite()){
    y= spheres_[rb.key()].center().y();
    } else {
    if (!rb.is_negative()) y= inf_;
    else y=-inf_;
    }*/
  return T::Point_2(t_.sphere(ra).center()[plane_coordinate(0).index()], 
		    t_.sphere(rb).center()[plane_coordinate(1).index()]);
}


/*Slice::T::Line_3 Slice::in_line_rz(Sds::Curve r, NT z) const {
  CGAL_precondition(r.is_rule());
  T::Point_3 pt(t_.center_c(r.key(),T::Coordinate_index(0)),
		t_.center_c(r.key(),T::Coordinate_index(1)),z);
  //switch (r.part() & (~Sds::Curve::IN_BIT)) {
  if (r.is_top())
    return T::Line_3(pt, T::Vector_3(0,-1,0));
  else if (r.is_bottom())
    return T::Line_3(pt, T::Vector_3(0,1,0));
  else if (r.is_left())
    return T::Line_3(pt, T::Vector_3(1,0,0));
  else if (r.is_right())
    return T::Line_3(pt, T::Vector_3(-1,0,0));
  else {
    CGAL_assertion(0);
    return T::Line_3();
  }
}
  
Slice::T::Line_3 Slice::out_line_rz(Sds::Curve r, NT z) const {
  CGAL_precondition(r.is_rule());
  T::Point_3 pt(t_.center_c(r.key(),T::Coordinate_index(0)),
		t_.center_c(r.key(),T::Coordinate_index(1)),z);
  //switch (r.part() & (~Sds::Curve::IN_BIT)) {
  if (r.is_top())
    return T::Line_3(pt, T::Vector_3(0,1,0));
  else if (r.is_bottom())
    return T::Line_3(pt, T::Vector_3(0,-1,0));
  else if (r.is_left())
    return T::Line_3(pt, T::Vector_3(-1,0,0));
  else if (r.is_right())
    return T::Line_3(pt, T::Vector_3(1,0,0));
  else {
    CGAL_assertion(0);
    return T::Line_3();
  }
  }*/

Slice::T::Line_3 Slice::positive_line_rz(T::Key k, Coordinate_index i, NT z) const {
  T::FT p[3];
  p[plane_coordinate(0).index()]= t_.sphere(k).center()[plane_coordinate(0).index()];
  p[plane_coordinate(1).index()]= t_.sphere(k).center()[plane_coordinate(1).index()];
  p[sweep_coordinate().index()]=z;
  T::Point_3 pt(p[0],p[1],p[2]);
  //switch (r.part() & (~Sds::Curve::IN_BIT)) {
  T::FT v[3]={0,0,0};
  if (i==plane_coordinate(0)) {
    v[plane_coordinate(1).index()]=1;
  } else {
    v[plane_coordinate(0).index()]=1;
  }
  return T::Line_3(pt, T::Vector_3(v[0], v[1], v[2]));
}
  
Slice::T::Line_3 Slice::negative_line_rz(T::Key k, Coordinate_index i, NT z) const {
  T::FT p[3];
  p[plane_coordinate(0).index()]= t_.sphere(k).center()[plane_coordinate(0).index()];
  p[plane_coordinate(1).index()]= t_.sphere(k).center()[plane_coordinate(1).index()];
  p[sweep_coordinate().index()]=z;
  T::Point_3 pt(p[0],p[1],p[2]);
  //switch (r.part() & (~Sds::Curve::IN_BIT)) {
  T::FT v[3]={0,0,0};
  if (i==plane_coordinate(0)) {
    v[plane_coordinate(1).index()]=-1;
  } else {
    v[plane_coordinate(0).index()]=-1;
  }
  return T::Line_3(pt, T::Vector_3(v[0], v[1], v[2]));
}



