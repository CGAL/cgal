#include <CGAL/Arrangement_of_spheres_3/Slice.h>

#define DPRINT(x) x





/* 
   2D -------------------------------------------------------------
*/

bool Slice::intersects_rz(T::Key a, NT z) const {
  CGAL::Comparison_result cmp= CGAL::compare(t_.sphere(a).squared_radius(),
					     CGAL::square(t_.center_c(a,T::Coordinate_index(2))-z));
  
  return cmp!= CGAL::SMALLER;
}


Slice::T::Circle_2 Slice::circle_rz(T::Key a, NT z) const {
  NT r2= t_.sphere(a).squared_radius()
    -  CGAL::square(t_.center_c(a,T::Coordinate_index(2))-z);
  CGAL_assertion(r2>=0);
  T::Circle_2  c(T::Point_2(t_.sphere(a).center().x(), 
			    t_.sphere(a).center().y()), r2);
  return c;
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
  if (pt.type()== Sds::Point::RR){
    T::Point_2 p= compute_rule_rule_intersection(pt.rule(0), pt.rule(1));
    T::Point_3 p3(p.x(), p.y(), 0);
    return  T::Sphere_point_3(p3, T::Line_3(p3,
					    T::Vector_3(0, 0, 1)));
  } else if (pt.type() == Sds::Point::SR) {
    if (pt.rule(0).key() == pt.sphere(0).key() || pt.rule(0).is_same_side(pt.sphere(0))) {
      T::Sphere_3 s(t_.sphere(pt.sphere(0).key()));
      T::Line_3 l(in_line(pt.rule(0), z));
      //std::cout << s << std::endl;
      //std::cout << l << std::endl;
      T::Sphere_point_3 sp(s, l);
      return sp;
    } else {
      return T::Sphere_point_3(t_.sphere(pt.sphere(0).key()), out_line(pt.rule(0), z));
    }
  } else if (pt.sphere(0).key() == pt.sphere(1).key()) {
    int ipt=static_cast<int>(pt.sphere(0).part() & pt.sphere(1).part()) 
      & (~Sds::Curve::ARC_BIT);
    Sds::Curve::Part cpt= static_cast<Sds::Curve::Part>(ipt);
    Sds::Curve rule(pt.sphere(0).key(), cpt);
    CGAL_assertion(rule.is_rule());
    return T::Sphere_point_3(t_.sphere(pt.sphere(0).key()), in_line(rule, z));
  } else {
    //std::cout << "Computing point for " << pt << std::endl;
    CGAL_precondition(pt.sphere(0).key() != pt.sphere(1).key());
    T::Point_2 cp= center_point_rz(pt.sphere(0).key(), pt.sphere(1).key(), z);
    T::Vector_3 v= t_.center(pt.sphere(0).key()) - t_.center(pt.sphere(1).key());
    //std::cout << "pt = " << cp << std::endl;
    //std::cout << "v = " << v << std::endl;
    T::Line_3 l(T::Point_3(cp.x(), cp.y(), z), T::Vector_3(-v.y(), v.x(), 0));
    
    T::Sphere_point_3 sli(t_.sphere(pt.sphere(1).key()), l);
    //std::cout << sli << std::endl;
    //std::cout << sli.approximate_coordinate(0) << ", " <<  sli.approximate_coordinate(1) << std::endl;
    T::Sphere_point_3 osli(t_.sphere(pt.sphere(1).key()), l);
    //std::cout << osli.approximate_coordinate(0) << ", " <<  osli.approximate_coordinate(1) << std::endl;
    CGAL_exactness_assertion(sli.compare(osli,Coordinate_index(0)) == CGAL::EQUAL);
    CGAL_exactness_assertion(sli.compare(osli,Coordinate_index(1)) == CGAL::EQUAL);
    CGAL_exactness_assertion(sli.compare(osli,Coordinate_index(2)) == CGAL::EQUAL);
    CGAL_assertion(sli.is_valid());
    CGAL_exactness_assertion(CGAL::abs(sli.exact_coordinate(Coordinate_index(0))) 
			     < t_.inf());
    CGAL_exactness_assertion(CGAL::abs(sli.exact_coordinate(Coordinate_index(1)))
			     < t_.inf());
    return sli;
  }
}

/* 
   constructions----------------------------------------------------
*/






Slice::T::Point_2 Slice::compute_rule_rule_intersection(Sds::Curve ra,
							Sds::Curve rb) const {
  CGAL_precondition(ra.is_rule() && rb.is_rule());
  CGAL_precondition(ra.is_vertical());
  CGAL_precondition(!rb.is_vertical());
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
  return T::Point_2(t_.center_c(ra.key(), T::Coordinate_index(0)), 
		    t_.center_c(rb.key(), T::Coordinate_index(1)));
}


Slice::T::Line_3 Slice::in_line(Sds::Curve r, NT z) const {
  CGAL_precondition(r.is_rule());
  T::Point_3 pt(t_.center_c(r.key(),T::Coordinate_index(0)),
		t_.center_c(r.key(),T::Coordinate_index(1)),z);
  switch (r.part() & (~Sds::Curve::IN_BIT)) {
  case Sds::Curve::T_RULE:
    return T::Line_3(pt, T::Vector_3(0,-1,0));
  case Sds::Curve::B_RULE:
    return T::Line_3(pt, T::Vector_3(0,1,0));
  case Sds::Curve::L_RULE:
    return T::Line_3(pt, T::Vector_3(1,0,0));
  case Sds::Curve::R_RULE:
    return T::Line_3(pt, T::Vector_3(-1,0,0));
  default:
    CGAL_assertion(0);
    return T::Line_3();
  }
}
  
Slice::T::Line_3 Slice::out_line(Sds::Curve r, NT z) const {
  CGAL_precondition(r.is_rule());
  T::Point_3 pt(t_.center_c(r.key(),T::Coordinate_index(0)),
		t_.center_c(r.key(),T::Coordinate_index(1)),z);
  switch (r.part() & (~Sds::Curve::IN_BIT)) {
  case Sds::Curve::T_RULE:
    return T::Line_3(pt, T::Vector_3(0,1,0));
  case Sds::Curve::B_RULE:
    return T::Line_3(pt, T::Vector_3(0,-1,0));
  case Sds::Curve::L_RULE:
    return T::Line_3(pt, T::Vector_3(-1,0,0));
  case Sds::Curve::R_RULE:
    return T::Line_3(pt, T::Vector_3(1,0,0));
  default:
    CGAL_assertion(0);
    return T::Line_3();
  }
}



