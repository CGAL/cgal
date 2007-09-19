#include <CGAL/Arrangement_of_spheres_3/Rational_cross_section.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
void
Rational_cross_section CGAL_AOS3_TARG ::audit() const {
  /*cs_.audit();
  for (CGAL_AOS3_TYPENAME CCS::Vertex_iterator vit= cs_.vertices_begin();
       vit != vs_.vertices_end(); ++vit) {
    if (vit->point().is_sphere_extremum()) {

    } else if (vit->
    }*/
  for (CGAL_AOS3_TYPENAME Traits::Sphere_3_key_const_iterator it= tr_.sphere_3_keys_begin();
       it != tr_.sphere_3_keys_end(); ++it){
    CGAL_AOS3_TYPENAME CCS::Halfedge_handle h= ccs_.a_halfedge(*it);
    if (h != CGAL_AOS3_TYPENAME CCS::Halfedge_handle()) {
      while (!h->vertex()->point().is_sphere_extremum() || h->vertex()->point().sphere_extremum_index() != Rule_direction(0)) {
	h= ccs_.next_edge_on_circle(h);
      }
      std::vector<CGAL_AOS3_TYPENAME CCS::Vertex_handle> tr, tl, bl, br;
      while (!h->vertex()->point().is_sphere_extremum() || h->vertex()->point().sphere_extremum_index() != Rule_direction(1)) {
	tr.push_back(h->vertex());
	h= ccs_.next_edge_on_circle(h);
      }
      while (!h->vertex()->point().is_sphere_extremum() || h->vertex()->point().sphere_extremum_index() != Rule_direction(2)) {
	tl.push_back(h->vertex());
	h= ccs_.next_edge_on_circle(h);
      }
      while (!h->vertex()->point().is_sphere_extremum() || h->vertex()->point().sphere_extremum_index() != Rule_direction(3)) {
	bl.push_back(h->vertex());
	h= ccs_.next_edge_on_circle(h);
      }
      while (!h->vertex()->point().is_sphere_extremum() || h->vertex()->point().sphere_extremum_index() != Rule_direction(0)) {
	br.push_back(h->vertex());
	h= ccs_.next_edge_on_circle(h);
      }
      std::vector<CGAL_AOS3_TYPENAME CCS::Vertex_handle> l[2][2];
      l[0][0].push_back(bl.front());
      l[0][0].insert(l[0][0].end(), tl.rbegin(), tl.rend());
      l[0][0].insert(l[0][0].end(), tr.rbegin(), tr.rend());
     
      l[0][1].insert(l[0][1].end(), bl.begin(), bl.end());
      l[0][1].insert(l[0][1].end(), br.begin(), br.end());
      l[0][1].push_back(tr.front());

      l[1][0].push_back(br.front());
      l[1][0].insert(l[1][0].end(), bl.rbegin(), bl.rend());
      l[1][0].insert(l[1][0].end(), tl.rbegin(), tl.rend());
     
      l[1][1].insert(l[1][1].end(), br.begin(), br.end());
      l[1][1].insert(l[1][1].end(), tr.begin(), tr.end());
      l[1][1].push_back(tl.front());

      for (int i=0; i< 2; ++i) {
	for (int j=0; j< 2; ++j) {
	  for (unsigned int k=0; k< l[i][j].size(); ++k){
	    CGAL_LOG(Log::LOTS,  l[i][j][k]->point() << " ");
	  }
	  CGAL_LOG(Log::LOTS, std::endl);

	  for (unsigned int k=1; k< l[i][j].size(); ++k){
	    CGAL_assertion(tr_.compare_c(sphere_point(l[i][j][k-1]->point()),
					 sphere_point(l[i][j][k]->point()),
					 plane_coordinate(i)) != CGAL::LARGER);
	  }
	}
      }
    }
  }

  for (CGAL_AOS3_TYPENAME CCS::Halfedge_const_iterator hit = ccs_.halfedges_begin();
       hit != ccs_.halfedges_end(); ++hit) {
    if (hit->curve().is_rule() && (hit->curve().is_inside() && hit->curve().is_vertical()
				   || !hit->curve().is_inside() && !hit->curve().is_vertical())) {
      CGAL_assertion(tr_.compare_c(sphere_point(hit->opposite()->vertex()->point()),
				   sphere_point(hit->vertex()->point()),
				   other_plane_coordinate(hit->curve().constant_coordinate())) != CGAL::LARGER);
    }
  }
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Rational_cross_section CGAL_AOS3_TARG ::Traits::Sphere_point_3 
Rational_cross_section CGAL_AOS3_TARG ::sphere_point(CGAL_AOS3_TYPENAME CCS::Point pt) const {
  if (pt.is_rule_rule()){
    CGAL_AOS3_TYPENAME Traits::Point_2 p= compute_rule_rule_intersection(pt.rule_key(plane_coordinate(0)),
									 pt.rule_key(plane_coordinate(1)));
    CGAL_AOS3_TYPENAME Traits::FT pt[3]={0,0,0};
    pt[plane_coordinate(0).index()]=p.x();
    pt[plane_coordinate(1).index()]=p.y();
    CGAL_AOS3_TYPENAME Traits::Point_3 p3(pt[0], pt[1], pt[2]);
  
    CGAL_AOS3_TYPENAME Traits::Sphere_point_3 ret(p3);
    if (!ret.is_valid()) {
      throw p3;
    } else return ret;
  } else if (pt.is_sphere_rule()) {
    CGAL_AOS3_TYPENAME Traits::Sphere_3 s(tr_.sphere_3(pt.sphere_key()));
    CGAL_AOS3_TYPENAME Traits::Line_3 l;
    if (pt.is_smaller()) {
      l= positive_line(pt.rule_key(), pt.rule_constant_coordinate());
    } else {
      l = negative_line(pt.rule_key(), pt.rule_constant_coordinate());
    }
    CGAL_AOS3_TYPENAME Traits::Sphere_point_3 sp (s,l);
    if (!sp.is_valid()) {
      std::cerr << "Constructing SR point for " << pt << " at " << z_ << std::endl;
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
      return CGAL_AOS3_TYPENAME Traits::Sphere_point_3(t_.sphere(pt.sphere(0).key()), 
      in_line_rz(rule, z));*/
    return CGAL_AOS3_TYPENAME Traits::Sphere_point_3();
  } else {
    //std::cout << "Computing point for " << pt << std::endl;
    CGAL_precondition(pt.sphere_key(0) != pt.sphere_key(1));
    CGAL_AOS3_TYPENAME Traits::Point_2 cp= center_point(pt.sphere_key(0), pt.sphere_key(1));
    CGAL_AOS3_TYPENAME Traits::Vector_3 v= tr_.sphere_3(pt.sphere_key(0)).center() 
      - tr_.sphere_3(pt.sphere_key(1)).center();
    //std::cout << "cp = " << cp << std::endl;
    //std::cout << "v = " << v << std::endl;
    CGAL_AOS3_TYPENAME Traits::FT p[3],vv[3];
    p[plane_coordinate(0).index()]= cp[0];
    p[plane_coordinate(1).index()]= cp[1];
    p[sweep_coordinate().index()]= z_;
    vv[plane_coordinate(0).index()]= -v[plane_coordinate(1).index()];
    vv[plane_coordinate(1).index()]= v[plane_coordinate(0).index()];
    vv[sweep_coordinate().index()]=0;
    //std::cout << p[0] << " " << p[1] << " " << p[2] << std::endl;
    //std::cout << vv[0] << " " << vv[1] << " " << vv[2] << std::endl;
    CGAL_AOS3_TYPENAME Traits::Line_3 l(CGAL_AOS3_TYPENAME Traits::Point_3(p[0], p[1], p[2]),
					CGAL_AOS3_TYPENAME Traits::Vector_3(vv[0], vv[1], vv[2]));
    
    CGAL_AOS3_TYPENAME Traits::Sphere_point_3 sli(tr_.sphere_3(pt.sphere_key(1)), l);
    if (!sli.is_valid()) {
      std::cerr << "Constructing point for " << pt << " at " << z_ << std::endl;
      std::cerr << "Line " << l << " does not intersect sphere " 
		<< tr_.sphere_3(pt.sphere_key(1)) << std::endl;
      throw l.projection(tr_.sphere_3(pt.sphere_key(1)).center());
    }
    //std::cout << sli << std::endl;
    //std::cout << sli.approximate_coordinate(0) << ", " <<  sli.approximate_coordinate(1) << std::endl;
    CGAL_AOS3_TYPENAME Traits::Sphere_point_3 osli(tr_.sphere_3(pt.sphere_key(1)), l);
    //std::cout << osli.approximate_coordinate(0) << ", " <<  osli.approximate_coordinate(1) << std::endl;
    
    CGAL_exactness_assertion(sli.compare(osli,Coordinate_index::X()) == CGAL::EQUAL);
    CGAL_exactness_assertion(sli.compare(osli,Coordinate_index::Y()) == CGAL::EQUAL);
    CGAL_exactness_assertion(sli.compare(osli,Coordinate_index::Z()) == CGAL::EQUAL);
    CGAL_assertion(sli.is_valid());
    CGAL_exactness_assertion(CGAL::abs(sli.exact_coordinate(plane_coordinate(0))) 
			     < tr_.max_coordinate());
    CGAL_exactness_assertion(CGAL::abs(sli.exact_coordinate(plane_coordinate(1)))
			     < tr_.max_coordinate());
    return sli;
  }
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Rational_cross_section CGAL_AOS3_TARG ::Traits::Point_2 
Rational_cross_section CGAL_AOS3_TARG ::rule_rule_intersection(CGAL_AOS3_TYPENAME Traits::Sphere_3_key ra,
		       CGAL_AOS3_TYPENAME Traits::Sphere_3_key rb) const {
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
  return CGAL_AOS3_TYPENAME Traits::Point_2(tr_.sphere_3(ra).center()[plane_coordinate(0).index()], 
					    tr_.sphere_3(rb).center()[plane_coordinate(1).index()]);
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Rational_cross_section CGAL_AOS3_TARG ::Traits::Circle_2
Rational_cross_section CGAL_AOS3_TARG ::circle(CGAL_AOS3_TYPENAME Traits::Sphere_3_key a) const {
  CGAL_AOS3_TYPENAME Traits::Sphere_3 s= tr_.sphere_3(a);
  CGAL_AOS3_TYPENAME Traits::FT r2= tr_.sphere_3(a).squared_radius()
    -  CGAL::square(tr_.sphere_3(a).center()[sweep_coordinate().index()]-z_);
  if (r2 >=0) {
    CGAL_AOS3_TYPENAME Traits::Circle_2  c(CGAL_AOS3_TYPENAME Traits::Point_2(tr_.sphere_3(a).center()[plane_coordinate(0).index()], 
									      tr_.sphere_3(a).center()[plane_coordinate(1).index()]), r2);
    return c;
  } else {
    return CGAL_AOS3_TYPENAME Traits::Circle_2(CGAL_AOS3_TYPENAME Traits::Point_2(tr_.sphere_3(a).center()[plane_coordinate(0).index()], 
										  tr_.sphere_3(a).center()[plane_coordinate(1).index()]), 0);
  }
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Rational_cross_section CGAL_AOS3_TARG ::Traits::Line_3 
Rational_cross_section CGAL_AOS3_TARG ::positive_line(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
						      Coordinate_index i) const {
  CGAL_AOS3_TYPENAME Traits::FT p[3];
  p[plane_coordinate(0).index()]= tr_.sphere_3(k).center()[plane_coordinate(0).index()];
  p[plane_coordinate(1).index()]= tr_.sphere_3(k).center()[plane_coordinate(1).index()];
  p[sweep_coordinate().index()]=z_;
  CGAL_AOS3_TYPENAME Traits::Point_3 pt(p[0],p[1],p[2]);
  //switch (r.part() & (~Sds::Curve::IN_BIT)) {
  CGAL_AOS3_TYPENAME Traits::FT v[3]={0,0,0};
  if (i==plane_coordinate(0)) {
    v[plane_coordinate(1).index()]=1;
  } else {
    v[plane_coordinate(0).index()]=1;
  }
  return CGAL_AOS3_TYPENAME Traits::Line_3(pt, CGAL_AOS3_TYPENAME Traits::Vector_3(v[0], v[1], v[2]));
}
  


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Rational_cross_section CGAL_AOS3_TARG ::Traits::Line_3 
Rational_cross_section CGAL_AOS3_TARG ::negative_line(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
						      Coordinate_index i) const {
  CGAL_AOS3_TYPENAME Traits::FT p[3];
  p[plane_coordinate(0).index()]= tr_.sphere_3(k).center()[plane_coordinate(0).index()];
  p[plane_coordinate(1).index()]= tr_.sphere_3(k).center()[plane_coordinate(1).index()];
  p[sweep_coordinate().index()]=z_;
  CGAL_AOS3_TYPENAME Traits::Point_3 pt(p[0],p[1],p[2]);
  //switch (r.part() & (~Sds::Curve::IN_BIT)) {
  CGAL_AOS3_TYPENAME Traits::FT v[3]={0,0,0};
  if (i==plane_coordinate(0)) {
    v[plane_coordinate(1).index()]=-1;
  } else {
    v[plane_coordinate(0).index()]=-1;
  }
  return CGAL_AOS3_TYPENAME Traits::Line_3(pt, CGAL_AOS3_TYPENAME Traits::Vector_3(v[0], v[1], v[2]));
}


CGAL_AOS3_TEMPLATE
bool Rational_cross_section CGAL_AOS3_TARG ::intersects(CGAL_AOS3_TYPENAME Traits::Sphere_3_key a) const {
  CGAL_AOS3_TYPENAME Traits::FT p[3];
  p[plane_coordinate(0).index()]=0;
  p[plane_coordinate(1).index()]=1;
  p[sweep_coordinate().index()]=z_;
  return tr_.bounded_side_of_sphere_c(CGAL_AOS3_TYPENAME Traits::Sphere_point_3(CGAL_AOS3_TYPENAME Traits::Point_3(p[0], p[1], p[2])),
				      a,
				      sweep_coordinate()); 
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Rational_cross_section CGAL_AOS3_TARG ::Traits::Point_2
Rational_cross_section CGAL_AOS3_TARG ::center_point(CGAL_AOS3_TYPENAME Traits::Sphere_3_key a, 
						     CGAL_AOS3_TYPENAME Traits::Sphere_3_key b) const {
  //CGAL_AOS3_TYPENAME Traits::Intersect_with_sweep is= tr_.intersect_with_sweep_object(z);
  CGAL_AOS3_TYPENAME Traits::Circle_2 ca= circle(a);
  CGAL_AOS3_TYPENAME Traits::Circle_2 cb= circle(b);
  CGAL_AOS3_TYPENAME Traits::Geom_traits::Compute_squared_length_2 csl= tr_.geom_traits_object().compute_squared_length_2_object();
  CGAL_AOS3_TYPENAME Traits::FT c02= csl(ca.center()-ORIGIN);
  CGAL_AOS3_TYPENAME Traits::FT c12= csl(cb.center()-ORIGIN);
  CGAL_AOS3_TYPENAME Traits::FT c0c1= (cb.center()-ORIGIN)*(ca.center()-ORIGIN);
  CGAL_AOS3_TYPENAME Traits::FT sqr=c02 + c12 -2*c0c1;
  CGAL_AOS3_TYPENAME Traits::FT t=(cb.squared_radius() - ca.squared_radius() +sqr)
    /(2*sqr);
  CGAL_AOS3_TYPENAME Traits::Point_2 p= ORIGIN+ (t*(ca.center()-ORIGIN) + (1-t)*(cb.center()-ORIGIN));
    
  return p;
}

CGAL_AOS3_END_INTERNAL_NAMESPACE
