#include <CGAL/Arrangement_of_spheres_3/Slice.h>

#define DPRINT(x) x







/* 
   2D -------------------------------------------------------------
*/

  Slice::T::Point_2 Slice::center_point_rz(int a, int b, NT z) const {
    T::Intersect_with_sweep is= tr_.intersect_with_sweep_object(z);
    T::Circle_2 ca= is(spheres_[a]);
    T::Circle_2 cb= is(spheres_[b]);
    T::Geometric_kernel::Compute_squared_length_2 csl= tr_.geometric_kernel_object().compute_squared_length_2_object();
    T::NT c02= csl(ca.center()-CGAL::ORIGIN);
    T::NT c12= csl(cb.center()-CGAL::ORIGIN);
    T::NT c0c1= (cb.center()-CGAL::ORIGIN)*(ca.center()-CGAL::ORIGIN);
    NT sqr=c02 + c12 -2*c0c1;
    NT t=(cb.squared_radius() - ca.squared_radius() +sqr)
      /(2*sqr);
    T::Point_2 p= CGAL::ORIGIN+ (t*(ca.center()-CGAL::ORIGIN) + (1-t)*(cb.center()-CGAL::ORIGIN));
    
    return p;
  }

  Slice::T::Point_2 Slice::display_point_rz(Sds::Point pt, NT z) const {
    T::Sphere_point_3 sp= sphere_point_rz(pt, z);
    //std::cout << "Exact point is " << sp << std::endl;
    return T::Point_2(sp.approximate_coordinate(0), sp.approximate_coordinate(1));
  }


Slice::T::Sphere_point_3 Slice::sphere_point_rz(Sds::Point pt, NT z) const {
  if (pt.type()== Sds::Point::RR){
    T::Point_2 p= compute_rule_rule_intersection(pt.rule(0), pt.rule(1));
    return  T::Sphere_point_3(p, T::Line_2(p, T::Vector_2(p.x(), 1)));
  } else if (pt.type() == Sds::Point::SR) {
    if (pt.rule(0).index() == pt.sphere(0).index() || pt.rule(0).is_same_side(pt.sphere(0))) {
      T::Sphere_3 s(spheres_[pt.sphere(0).index()]);
      T::Line_3 l(in_line(pt.rule(0), z));
      //std::cout << s << std::endl;
      //std::cout << l << std::endl;
      T::Sphere_point_3 sp(s, l);
      return sp;
    } else {
      return T::Sphere_point_3(spheres_[pt.sphere(0).index()], out_line(pt.rule(0), z));
    }
  } else if (pt.sphere(0).index() == pt.sphere(1).index()) {
    int ipt=static_cast<int>(pt.sphere(0).part() & pt.sphere(1).part()) 
      & (~Sds::Curve::ARC_BIT);
    Sds::Curve::Part cpt= static_cast<Sds::Curve::Part>(ipt);
    Sds::Curve rule(pt.sphere(0).index(), cpt);
    CGAL_assertion(rule.is_rule());
    return T::Sphere_point_3(spheres_[pt.sphere(0).index()], in_line(rule, z));
  } else {
    //std::cout << "Computing point for " << pt << std::endl;
    CGAL_precondition(pt.sphere(0).index() != pt.sphere(1).index());
    T::Point_2 cp= center_point_rz(pt.sphere(0).index(), pt.sphere(1).index(), z);
    T::Vector_3 v= spheres_[pt.sphere(0).index()].center() - spheres_[pt.sphere(1).index()].center();
    //std::cout << "pt = " << cp << std::endl;
    //std::cout << "v = " << v << std::endl;
    T::Line_3 l(T::Point_3(cp.x(), cp.y(), z), T::Vector_3(-v.y(), v.x(), 0));
    
    T::Sphere_point_3 sli(spheres_[pt.sphere(1).index()], l);
    //std::cout << sli << std::endl;
    //std::cout << sli.approximate_coordinate(0) << ", " <<  sli.approximate_coordinate(1) << std::endl;
    T::Sphere_point_3 osli(spheres_[pt.sphere(1).index()], l);
    //std::cout << osli.approximate_coordinate(0) << ", " <<  osli.approximate_coordinate(1) << std::endl;
    CGAL_exactness_assertion(sli.compare(osli,0) == CGAL::EQUAL);
    CGAL_exactness_assertion(sli.compare(osli,1) == CGAL::EQUAL);
    CGAL_exactness_assertion(sli.compare(osli,2) == CGAL::EQUAL);
    CGAL_assertion(sli.is_valid());
    CGAL_exactness_assertion(CGAL::abs(sli.exact_coordinate(0)) < inf_);
    CGAL_exactness_assertion(CGAL::abs(sli.exact_coordinate(1)) < inf_);
    return sli;
  }
}

/* 
   constructions----------------------------------------------------
*/



Slice::T::Plane_3 Slice::separating_plane(int a, int b) const {
  T::Vector_3 d(spheres_[b].center()-spheres_[a].center());
  T::Vector_3 n(-d[1], d[0], 0);
  return T::Plane_3(spheres_[b].center(), n);
}



Slice::T::Plane_3 Slice::equipower_plane(int a, int b) const {
  CGAL_precondition(a != b);
  T::Vector_3 n=2*(spheres_[a].center()-spheres_[b].center());
  NT d= disc(b) - disc(a);
  return T::Plane_3(n[0], n[1], n[2], d);
}

Slice::NT Slice::disc(int i) const {
  T::Vector_3 v= spheres_[i].center()-CGAL::ORIGIN;
  return v*v - spheres_[i].squared_radius();
}
  

Slice::NT Slice::center_c(int sphere, int C) const {
  //CGAL_precondition(c.is_rule());
  if  (sphere>=0) {
    return spheres_[sphere].center()[C];
  } else if (sphere==-1) return -inf_;
  else return inf_;
}


Slice::T::Sphere_3 sphere(int ind) const {
  return spheres_[ind];
}

Slice::T::Event_point_3 Slice::sphere_start(int ind) const {
  return T::Event_point_3 ep(spheres_[ind], T::Line(center(ind), T::Vector(0,0,1)));
}


Slice::T::Point_2 Slice::compute_rule_rule_intersection(Sds::Curve ra,
							Sds::Curve rb) const {
  CGAL_precondition(ra.is_rule() && rb.is_rule());
  CGAL_precondition(ra.is_vertical());
  CGAL_precondition(!rb.is_vertical());
  /*NT x,y;
  if (ra.is_finite()){
    x= spheres_[ra.index()].center().x();
  } else {
    if (!ra.is_negative()) x= inf_;
    else x=-inf_;
  }
  if (rb.is_finite()){
    y= spheres_[rb.index()].center().y();
  } else {
    if (!rb.is_negative()) y= inf_;
    else y=-inf_;
    }*/
  return T::Point_2(center_c(ra.index()), center_c(rb.index()));
}


Slice::T::Line_3 Slice::in_line(Sds::Curve r, NT z) const {
  CGAL_precondition(r.is_rule());
  T::Point_3 pt(spheres_[r.index()].center().x(),spheres_[r.index()].center().y(),z);
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
  T::Point_3 pt(spheres_[r.index()].center().x(),spheres_[r.index()].center().y(),z);
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



