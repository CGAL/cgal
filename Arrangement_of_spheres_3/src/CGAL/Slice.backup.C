#include <CGAL/Arrangement_of_spheres_3/Slice.h>
#include <CGAL/IO/Qt_examiner_viewer_2.h>

#define DPRINT(x) x


/*template <class It>
  Slice::Vertex_const_handle Slice::find_vertex(It b, It e) const {
  Vertex_const_handle v0=(*b)->vertex(), v1= (*b)->opposite()->vertex();
  for (; b != e; ++b){
  Vertex_const_handle v= (*b)->vertex();
  Vertex_const_handle ov= (*b)->opposite()->vertex();
  CGAL_assertion(v0== v || v1 == v
  || v0 == ov || v1 == ov);
  if (v0== v && v1 != ov || v0 == ov && v1 != v) return v0;
  else if (v1== v && v0 != ov || v1 == ov && v0 != v) return v1;
  }
  CGAL_assertion(0);
  return Vertex_const_handle();
  }*/

void Slice::draw_rz(Qt_examiner_viewer_2 *qtv, NT z) {
  *qtv << CGAL::RED;
  qtv->set_updating_box(true);
  T::Intersect_with_sweep is=tr_.intersect_with_sweep_object(z);
  /*std::vector<bool> active(spheres_.size(), false);
    for (Slice_data_structure::Halfedge_const_iterator hit= sds_.halfedges_begin();
    hit != sds_.halfedges_end(); ++hit){
    if (hit->curve().is_finite()) {
    active[hit->curve().index()]=true;
    }
    }
    *qtv << CGAL::BLACK;
   
    for (unsigned int i=0; i< active.size(); ++i){
    if (active[i]){
    *qtv << is(spheres_[i]);
    }
    }*/
    
    
  for (Slice_data_structure::Halfedge_const_iterator hit= sds_.halfedges_begin();
       hit != sds_.halfedges_end(); ++hit){
    bool marked= (marked_faces_.find(hit->face()) != marked_faces_.end() 
		  || marked_faces_.find(hit->opposite()->face()) != marked_faces_.end())
      || (marked_edges_.find(hit) != marked_edges_.end()
	  || marked_edges_.find(hit->opposite()) != marked_edges_.end());

    if (hit->curve().is_rule() && hit->curve().is_inside()){
      qtv->set_updating_box(false);
      //std::cout << "Displaying rule " << hit->curve() << std::endl;
      T::Point_2 t= display_point_rz(hit->vertex()->point(), z);
      T::Point_2 s= display_point_rz(hit->opposite()->vertex()->point(), z);
      if (marked) {
	*qtv << CGAL::RED;
      } else {
	*qtv << CGAL::GRAY;
      }
      *qtv << T::Geometric_kernel::Segment_2(t,s);
    } else if (hit->curve().is_arc() && hit->curve().is_inside()){
      qtv->set_updating_box(true);
      //std::cout << "Displaying arc " << hit->curve() << std::endl;
      T::Point_2 t= display_point_rz(hit->vertex()->point(), z);
      T::Point_2 s= display_point_rz(hit->opposite()->vertex()->point(), z);
      T::Circle_2 c= is(spheres_[hit->curve().index()]);
      if (marked) {
	*qtv << CGAL::RED;
      } else {
	*qtv << CGAL::BLACK;
      }
      qtv->new_circular_arc(c, s, t);
    }
  }
   
  for (Slice_data_structure::Vertex_const_iterator hit= sds_.vertices_begin();
       hit != sds_.vertices_end(); ++hit){
    bool marked= (marked_vertices_.find(hit) != marked_vertices_.end());
    T::Point_2 p= display_point_rz(hit->point(), z);
    if (marked) {
      *qtv << CGAL::RED;
    } else {
      *qtv << CGAL::BLUE;
    }
    if (hit->point().is_finite()) {
      qtv->set_updating_box(true);
    } else {
      qtv->set_updating_box(false);
    }
    *qtv << p;
    //int sa= hit->point().first().index();
    //int sb= hit->point().second().index();
    std::ostringstream out;
    //if (sa != sb) {
    out << hit->point().first() << ":" << hit->point().second();
    //} else {
    //	out << sa;
    //}
    *qtv << CGAL::GRAY;
    *qtv << out.str().c_str();
    //    qtv->redraw();
  }
}









void Slice::set_rz(NT z) {
  Slice_arrangement sa(spheres_.begin(),spheres_.end(), z, inf_);
  sds_.clear();
  sds_.set_is_building(true);

  for (Slice_arrangement::Face_iterator fit= sa.faces_begin(); fit != sa.faces_end(); ++fit){
    sds_.new_face(fit->begin(), fit->end());
  }
    
  sds_.set_is_building(false);
}












void Slice::point_sphere_orientation(const T::Sphere_location &sl,
			     int sphere,
			     std::vector<int> &locations
			     /*,
			       std::vector<Sds::Curve> &edges*/) const {
  typedef T::Sphere_location SL;
  if (locations[sphere]==0){
    locations[sphere]=sl(spheres_[sphere]);
    DPRINT(std::cout << "For sphere " << sphere << " got " << SL::decode(locations[sphere]) << std::endl);
    /*int l= locations[sphere];
      if ((l & SL::IN_BIT) && (l & SL::OUT_BIT)) {
      if (l & SL::B_BIT) {
      if (l & SL::L_BIT) {
      edges.push_back(Sds::Curve(sphere, Sds::Curve::LB_ARC));
      } else {
      CGAL_assertion(l &SL::R_BIT);
      edges.push_back(Sds::Curve(sphere, Sds::Curve::RB_ARC));
      }
      } else {
      CGAL_assertion( l& SL::T_BIT);
      if (l & SL::L_BIT) {
      edges.push_back(Sds::Curve(sphere, Sds::Curve::LT_ARC));
      } else {
      CGAL_assertion(l &SL::R_BIT);
      edges.push_back(Sds::Curve(sphere, Sds::Curve::LT_ARC));
      }
      }
      } else if (l &SL::OUT_BIT) {
      if ((l & SL::T_BIT) && (l &SL::B_BIT) && (l& SL::L_BIT)){
      edges.push_back(Sds::Curve(sphere, Sds::Curve::L_RULE));
      } else if ((l & SL::T_BIT) && (l &SL::B_BIT) && (l& SL::R_BIT)){
      edges.push_back(Sds::Curve(sphere, Sds::Curve::R_RULE));
      } else if ((l & SL::L_BIT) && (l &SL::R_BIT) && (l& SL::T_BIT)){
      edges.push_back(Sds::Curve(sphere, Sds::Curve::T_RULE));
      } else if ((l & SL::L_BIT) && (l &SL::R_BIT) && (l& SL::B_BIT)){
      edges.push_back(Sds::Curve(sphere, Sds::Curve::L_RULE));
      }
      }*/
  }
}










bool Slice::locate_point_check_face(Face_const_iterator it,  const T::Sphere_location &sl,
				    std::vector<int> &locations) const {
  Halfedge_const_handle h= it->halfedge();
  bool finite=false;
  do {
    if (h->curve().is_finite()) {
      int sphere= h->curve().index();
      point_sphere_orientation(sl, sphere, locations);
      if (!h->curve().is_compatible_location(locations[sphere])) {
	//DPRINT(std::cout << "Nixed by edge " << h->curve() << std::endl);
	return false;
      }
      finite=true;
    }
    h= h->next();
  } while (h != it->halfedge());
  /*if (finite) {
    std::cout << "Face " << std::endl;;
    do {
    std::cout << h->curve() << std::endl;
    h= h->next();
    } while (h != it->halfedge());
    std::cout << std::endl;
    }*/

  return finite;
}










bool Slice::locate_point_check_face_vertices(T::Event_point_3 ep,
					     Face_const_iterator it) const {
  Halfedge_const_handle h= it->halfedge();
  do {
    if (h->vertex()->point().type() == Sds::Point::SS 
	&& h->curve().index() != h->next()->curve().index()) {
      // to handle tangencies use the actual curves, this gets the
      //orientation correct
      //Sds::Point npt(h->next()->curve(), h->curve());
      Sds::Point npt= h->vertex()->point();
      /*if (npt != npt){
	std::cout << "Rearranging vertex from " << h->vertex()->point() 
		  << " to " << npt << std::endl;
		  }*/
      if (oriented_side_of_center_plane(npt.sphere(0).index(),
					npt.sphere(1).index(),
					T::Point_2(ep.simple_coordinate(0),
						   ep.simple_coordinate(1)))
	  == CGAL::ON_NEGATIVE_SIDE) {
	std::cout << "Face nixed by vertex " << npt << std::endl;
	return false;
      }
    }
    h= h->next();
  } while (h != it->halfedge());

  return true;
}









 

Slice::Face_const_handle Slice::locate_point(T::Event_point_3 ep) const {
  if (CGAL::abs(ep.simple_coordinate(0)) > inf_
      || CGAL::abs(ep.simple_coordinate(1)) > inf_){
    std::cerr << "Coordinate out of range." << std::endl;
    CGAL_assertion(0);
  }
  std::vector<int> locations(spheres_.size(), 0);
  std::vector<Face_const_handle> faces;
  std::vector<Sds::Curve> edges;
  T::Sphere_location sl= tr_.sphere_location_object(ep);
  for (Face_const_iterator fit = sds_.faces_begin(); fit != sds_.faces_end(); ++fit){
    {
      std::cout << "Trying face ";
      Halfedge_const_handle h= fit->halfedge();
      do {
	std::cout << h->curve() << "--" << h->vertex()->point() << "--";
	h= h->next();
      } while (h != fit->halfedge());
      std::cout << std::endl;
    }
    bool ok=locate_point_check_face(fit, sl, locations/*, edges*/);
    if (ok) faces.push_back(fit);
  }

  if (faces.size()==1) {
    return  faces[0];
  } else if (faces.size() > 1) {
    //std::cout << "simplify this " << std::endl;
    std::vector<Face_const_handle> clean_faces;
    for (unsigned int i=0; i< faces.size(); ++i){
      if (true) {
	clean_faces.push_back(faces[i]);
	{
	  std::cout << "Faces cleaned " << std::endl;
	  Halfedge_const_handle h= faces[i]->halfedge();
	  do {
	    std::cout << h->curve() << "--";
	    h= h->next();
	  } while (h != faces[i]->halfedge());
	  std::cout << std::endl;
	}
      } else {
	std::cout << "Face unclean " << std::endl;
	Halfedge_const_handle h= faces[i]->halfedge();
	do {
	  std::cout << h->curve() << "--";
	  h= h->next();
	} while (h != faces[i]->halfedge());
	std::cout << std::endl;
      }
    }
    std::swap(faces, clean_faces);
     
  }

  if (faces.size() ==1) {
    return faces[0];
  } else if (faces.size() ==2) {
    Halfedge_const_handle h= faces[0]->halfedge();
    do {
      if (h->opposite()->face() == faces[1]) throw On_edge_exception(h);
      h= h->next();
    } while (h != faces[0]->halfedge());

    if (locate_point_check_face_vertices(ep, faces[0])) {
      std::cout << "Using center line to pick face." << std::endl;
      return faces[0];
    } else {
      CGAL_assertion( locate_point_check_face_vertices(ep, faces[1]));
      std::cout << "Using center line to pick face." << std::endl;
      return faces[1];
    }
    // 
    {
      std::cout << "Face 0 ";
      Halfedge_const_handle h= faces[0]->halfedge();
      do {
	std::cout << h->curve() << "--" << h->vertex()->point() << "--";
	h= h->next();
      } while (h != faces[0]->halfedge());
      std::cout << std::endl;
    }
    {
      std::cout << "\n\nFace 1 ";
      Halfedge_const_handle h= faces[1]->halfedge();
      do {
	std::cout << h->curve() << "--" << h->vertex()->point() << "--";
	h= h->next();
      } while (h != faces[1]->halfedge());
      std::cout << std::endl;
    }
  } else {
    CGAL_assertion(!faces.empty());
    Halfedge_const_handle h= faces[0]->halfedge();
    do {
      bool ok=true;
      for (unsigned int i=1; i< faces.size(); ++i) {
	if (sds_.has_vertex(faces[i], h->vertex())) {

	} else { ok=false;}
      }
      if (ok)   throw On_vertex_exception(h->vertex());
      h= h->next();
    } while (h != faces[0]->halfedge());
  }
   
  CGAL_assertion(0);
  return Face_const_handle();
}














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
      std::cout << s << std::endl;
      std::cout << l << std::endl;
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
  

  Slice::NT Slice::constant_coordinate(Sds::Curve c) const {
    CGAL_precondition(c.is_rule());
    if  (c.is_finite()) {
      return spheres_[c.index()].center()[c.constant_coordinate()];
    } else if (c.is_negative()) return -inf_;
    else return inf_;
  }

  Slice::T::Point_2 Slice::compute_rule_rule_intersection(Sds::Curve ra,
							  Sds::Curve rb) const {
    CGAL_precondition(ra.is_rule() && rb.is_rule());
    CGAL_precondition(ra.is_vertical());
    CGAL_precondition(!rb.is_vertical());
    NT x,y;
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
    }
    return T::Point_2(x,y);
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



/* 
   predictes--------------------------------------------------------
*/

CGAL::Comparison_result Slice::compare_sphere_to_plane(int s, NT c, int C) const {
  NT d= spheres_[s].center()[C]-c;
  if (CGAL::square(d) < spheres_[s].squared_radius()) return CGAL::EQUAL;
  else return compare_sphere_center_to_plane(s, c, C);
}

CGAL::Comparison_result  Slice::compare_sphere_center_to_plane(int s, NT c, int C) const {
  return CGAL::compare(spheres_[s].center()[C], c);
}

CGAL::Comparison_result  Slice::compare_equipower_point_to_plane(int a, int b, NT coord, int C) const {
  // can improve
  T::Plane_3 eqp= equipower_plane(a,b);
  T::Line_3 l(spheres_[a].center(), (spheres_[a].center()-spheres_[b].center()));
  CGAL::Object o= CGAL::intersection(eqp, l);
  T::Point_3 pt;
  CGAL_assertion(CGAL::assign(pt, o));
  CGAL::assign(pt, o);
  return CGAL::compare(pt[C], coord);
}

CGAL::Oriented_side Slice::oriented_side_of_center_plane(int a, int b,
							 T::Point_2 pt) const {
  T::Vector_3 d(spheres_[b].center()-spheres_[a].center());
  T::Line_2 l(T::Point_2(spheres_[a].center().x(),
			 spheres_[a].center().y()), 
	      T::Vector_2(d.x(), d.y()));
  return l.oriented_side(pt);
}
