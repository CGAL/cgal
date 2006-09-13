#include <CGAL/Arrangement_of_spheres_3/Slice.h>
#include <CGAL/IO/Qt_examiner_viewer_2.h>

#define DPRINT(x) x

std::ostream &Slice::write(Vertex_const_handle v, std::ostream &out) const {
  out << v->point();
  return out;
}
std::ostream &Slice::write(Halfedge_const_handle v, std::ostream &out) const {
  out << v->curve();
  return out;
}
std::ostream &Slice::write(Face_const_handle v, std::ostream &out) const {
  Halfedge_const_handle h=v->halfedge();
  do {
    out << h->curve() << "--" << h->vertex()->point() << "--";
    h=h->next();
  } while (h != v->halfedge());
  return out;
}


Slice::T::Key Slice::debug_new_sphere(T::Sphere_3 s) {
  T::Key r= t_.new_sphere(s);
  sds_.debug_add_sphere();
  //rule_events_.resize(rule_events_.size()+1);
  //CGAL_postcondition(halfedges_.size()== t_.number_of_spheres());

  return r;
}


Slice::DT::Point_2 Slice::display_point_rz(Sds::Point pt, NT z) const {
  try {
    T::Sphere_point_3 sp= sphere_point_rz(pt, z);
    //std::cout << "Exact point is " << sp << std::endl;
    return DT::Point_2(sp.approximate_coordinate(plane_coordinate(0)),
		       sp.approximate_coordinate(plane_coordinate(1)));
  } catch (T::Point_3 pt) {
    std::cout << "Point " << pt << " is no longer valid at " << z << std::endl;
    return DT::Point_2(pt[plane_coordinate(0).index()],
		       pt[plane_coordinate(1).index()]);
  }
}





void Slice::draw_rz(Qt_examiner_viewer_2 *qtv, NT z) {
  //t_.set_temp_sphere(T::Sphere_3(T::Point_3(0,0,z), 0));

  *qtv << CGAL::RED;
  qtv->set_updating_box(true);
  //T::Intersect_with_sweep is=t_.sphere_intersects_rule(z);
    
    
  /*for (T::Sphere_key_iterator sit= t_.sphere_keys_begin(); 
       sit != t_.sphere_keys_end(); ++sit){
    if (intersects_rz(*sit, z) {
      T::Circle_2 c2= circle_rz(*sit, z);
     
	c2= T::Circle_2(c2.center(), c2.squared_radius()*NT(1.01));
      } else {
	c2= T::Circle_2(c2.center(), c2.squared_radius()*NT(.99));
      }
      if (t_.sphere(*sit).center().z() != z){
	*qtv << CGAL::YELLOW;
	*qtv << c2;
      }
    }
    }*/

  for (Slice_data_structure::Halfedge_const_iterator hit= sds_.halfedges_begin();
       hit != sds_.halfedges_end(); ++hit){
    if (hit->curve().key().is_target()) continue;
    if (hit->curve().is_rule() && hit->curve().is_inside()){
      qtv->set_updating_box(false);
      //std::cout << "Displaying rule " << hit->curve() << std::endl;
      DT::Point_2 t= display_point_rz(hit->vertex()->point(), z);
      DT::Point_2 s= display_point_rz(hit->opposite()->vertex()->point(), z);
   
      *qtv << CGAL::GRAY;
      *qtv << DT::Segment_2(t,s);
    } else if (hit->curve().is_arc() && hit->curve().is_inside()){
      qtv->set_updating_box(true);
      //std::cout << "Displaying arc " << hit->curve() << std::endl;
      DT::Point_2 t= display_point_rz(hit->vertex()->point(), z);
      DT::Point_2 s= display_point_rz(hit->opposite()->vertex()->point(), z);
      //DT::Circle_2 c= ;
      if (t_.compare_sphere_center_c(hit->curve().key(), z,
				     sweep_coordinate())== CGAL::LARGER) {
	*qtv << CGAL::Color(150,50,50);
      } else {
	*qtv << CGAL::Color(50,150,50);
      }

      T::Circle_2 c2= circle_rz(hit->curve().key(), z);
      qtv->new_circular_arc(c2, s, t);
      
    }
  }
   
  for (Slice_data_structure::Vertex_const_iterator hit= sds_.vertices_begin();
       hit != sds_.vertices_end(); ++hit){
    if (!sds_.is_in_slice(hit)) continue;
    DT::Point_2 p= display_point_rz(hit->point(), z);
      *qtv << CGAL::BLUE;
    if (hit->point().is_finite()) {
      qtv->set_updating_box(true);
    } else {
      qtv->set_updating_box(false);
    }
    *qtv << p;
   
    std::ostringstream out;
    out << hit->point();
    /*if (hit->point().first().key() == hit->point().second().key()){
      out << hit->point().rule(0);
    } else {
      if (hit->point().first().is_arc()){
	out << hit->point().first().key();
      } else {
	out << hit->point().first();
      }
      out << ":";
      if (hit->point().second().is_arc()){
	out << hit->point().second().key();
      } else {
	out << hit->point().second();
      }
      }*/
    //out << hit->point().first() << ":" << hit->point().second();
    
    *qtv << CGAL::GRAY;
    *qtv << out.str().c_str();
  }

  *qtv << CGAL::Color(255, 155, 155);
  for (Intersections_2::const_iterator it= intersections_2_.begin();
       it != intersections_2_.end(); ++it) {
    if (it->second.first != Event_key()) {
      T::Event_point_3 ep= sim_->event_time(it->second.first);
      DT::Point_2 dp(ep.approximate_coordinate(plane_coordinate(0)), 
		     ep.approximate_coordinate(plane_coordinate(1)));
      *qtv << dp;
    }
    if (it->second.second != Event_key()) {
      T::Event_point_3 ep= sim_->event_time(it->second.second);
      DT::Point_2 dp(ep.approximate_coordinate(plane_coordinate(0)), 
		     ep.approximate_coordinate(plane_coordinate(1)));
      *qtv << dp;
    }
  }

  *qtv << CGAL::Color(255, 255, 255);
  for (Intersections_3::const_iterator it= intersections_3_.begin();
       it != intersections_3_.end(); ++it) {
    if (it->second.second != Event_key()) {
      {
	T::Event_point_3 ep= sim_->event_time(it->second.second);
	DT::Point_2 dp(ep.approximate_coordinate(plane_coordinate(0)), 
		       ep.approximate_coordinate(plane_coordinate(1)));
	*qtv << dp;
      }
      {
	T::Event_point_3 ep= sim_->event_time(it->second.second);
	DT::Point_2 dp(ep.approximate_coordinate(plane_coordinate(0)), 
		       ep.approximate_coordinate(plane_coordinate(1)));
	*qtv << dp;
      }
    }
  }
}


void Slice::draw_rz(CGAL::Qt_widget *qtv, NT z) {
  //t_.set_temp_sphere(T::Sphere_3(T::Point_3(0,0,z), 0));

  *qtv << CGAL::RED;
  //T::Intersect_with_sweep is=t_.sphere_intersects_rule(z);
    
    
  /*for (T::Sphere_key_iterator sit= t_.sphere_keys_begin(); 
       sit != t_.sphere_keys_end(); ++sit){
    if (intersects_rz(*sit, z) {
      T::Circle_2 c2= circle_rz(*sit, z);
     
	c2= T::Circle_2(c2.center(), c2.squared_radius()*NT(1.01));
      } else {
	c2= T::Circle_2(c2.center(), c2.squared_radius()*NT(.99));
      }
      if (t_.sphere(*sit).center().z() != z){
	*qtv << CGAL::YELLOW;
	*qtv << c2;
      }
    }
    }*/

  for (Slice_data_structure::Halfedge_const_iterator hit= sds_.halfedges_begin();
       hit != sds_.halfedges_end(); ++hit){
    if (hit->curve().key().is_target()) continue;
    if (hit->curve().is_rule() && hit->curve().is_inside()){
      //std::cout << "Displaying rule " << hit->curve() << std::endl;
      DT::Point_2 t= display_point_rz(hit->vertex()->point(), z);
      DT::Point_2 s= display_point_rz(hit->opposite()->vertex()->point(), z);
      if (hit->event() != Simulator::Event_key()) {
	if (hit->event() == sim_->null_event()) {
	  *qtv << CGAL::Color(150,50,50);
	} else {
	  *qtv << CGAL::Color(250,50,50);
	}
      } else {
	*qtv << CGAL::Color(200,200,200);
      }
      *qtv << DT::Segment_2(t,s);
    } else if (hit->curve().is_arc() && hit->curve().is_inside()){
      //std::cout << "Displaying arc " << hit->curve() << std::endl;
      DT::Point_2 t= display_point_rz(hit->vertex()->point(), z);
      DT::Point_2 s= display_point_rz(hit->opposite()->vertex()->point(), z);
      //DT::Circle_2 c= ;
      if (hit->event() != Simulator::Event_key()) {
	if (hit->event() == sim_->null_event()) {
	  *qtv << CGAL::Color(150,50,50);
	} else {
	  *qtv << CGAL::Color(250,50,50);
	}
      } else {
	*qtv << CGAL::Color(50,50,50);
      }
      
      T::Circle_2 c2= circle_rz(hit->curve().key(), z);
      *qtv << c2;
    }
  }
   
  for (Slice_data_structure::Vertex_const_iterator hit= sds_.vertices_begin();
       hit != sds_.vertices_end(); ++hit){
    if (!sds_.is_in_slice(hit)) continue;
    DT::Point_2 p= display_point_rz(hit->point(), z);
      *qtv << CGAL::BLUE;
 
    *qtv << p;
   
    std::ostringstream out;
    out << hit->point();
    
  }

  /*

   if (numv < 150) {
      std::ostringstream oss;
      oss << *it;
      w->get_painter().drawText(w->x_pixel(CGAL::to_double(pt.x()))+3,
				w->y_pixel(CGAL::to_double(pt.y()))-3,
				QString(oss.str().c_str()));
    }
  */

  //std::set<Intersection_3> certificates_3_;

  draw_events_rz(qtv, z);
}



void Slice::draw_marked_rz(Qt_examiner_viewer_2 *qtv, NT z) {
  //t_.set_temp_sphere(T::Sphere_3(T::Point_3(0,0,z), 0));

  *qtv << CGAL::RED;
  qtv->set_updating_box(true);
  //T::Intersect_with_sweep is=t_.sphere_intersects_rule(z);
    
    
  for (Slice_data_structure::Halfedge_const_iterator hit= sds_.halfedges_begin();
       hit != sds_.halfedges_end(); ++hit){
    bool marked= (marked_faces_.find(hit->face()) != marked_faces_.end() 
		  || marked_faces_.find(hit->opposite()->face()) != marked_faces_.end())
      || (marked_edges_.find(hit) != marked_edges_.end()
	  || marked_edges_.find(hit->opposite()) != marked_edges_.end());
    if (!marked) continue;

    if (hit->curve().is_rule() && hit->curve().is_inside()){
      qtv->set_updating_box(false);
      //std::cout << "Displaying rule " << hit->curve() << std::endl;
      DT::Point_2 t= display_point_rz(hit->vertex()->point(), z);
      DT::Point_2 s= display_point_rz(hit->opposite()->vertex()->point(), z);
     
      *qtv << CGAL::RED;
     
      *qtv << DT::Segment_2(t,s);
    } else if (hit->curve().is_arc() && hit->curve().is_inside()){
      qtv->set_updating_box(true);
      //std::cout << "Displaying arc " << hit->curve() << std::endl;
      DT::Point_2 t= display_point_rz(hit->vertex()->point(), z);
      DT::Point_2 s= display_point_rz(hit->opposite()->vertex()->point(), z);
      //DT::Circle_2 c= ;
     
      *qtv << CGAL::RED;
      
      qtv->new_circular_arc(circle_rz(hit->curve().key(), z), s, t);
    }
  }
   
  for (Slice_data_structure::Vertex_const_iterator hit= sds_.vertices_begin();
       hit != sds_.vertices_end(); ++hit){
    bool marked= (marked_vertices_.find(hit) != marked_vertices_.end());
    if (!marked) continue;

    DT::Point_2 p= display_point_rz(hit->point(), z);
  
      *qtv << CGAL::RED;
   
    if (hit->point().is_finite()) {
      qtv->set_updating_box(true);
    } else {
      qtv->set_updating_box(false);
    }
    *qtv << p;
   
    std::ostringstream out;
    Sds::Point pt= hit->point();
    
    if (pt.is_sphere_extremum()) {
      out << pt.sphere_key();
    } else if (pt.is_sphere_rule()) {
      out << pt.sphere_key() << ":" << pt.rule_key();
    } else if (pt.is_rule_rule()) {
      out << pt.rule_key(plane_coordinate(0)) << ":" << pt.rule_key(plane_coordinate(1));
    } else {
      out << pt.sphere_key(0) << ":" << pt.sphere_key(1);
    }
    /*if (hit->point().first().is_arc()){
      out << hit->point().first().key();
    } else {
      out << hit->point().first();
    }
    out << ":";
    if (hit->point().second().is_arc()){
      out << hit->point().second().key();
    } else {
      out << hit->point().first();
      }*/
    //out << hit->point().first() << ":" << hit->point().second();
    
    *qtv << CGAL::GRAY;
    *qtv << out.str().c_str();
  }
}









void Slice::initialize_at(NT z) {
  sim_->clear();
  sim_->set_current_time(Simulator::Time(z));
  Slice_arrangement sa(t_.spheres_begin(),t_.spheres_end(), z, 
		       t_.max_coordinate());
  sds_.clear();
  sds_.set_is_building(true);

  for (Slice_arrangement::Face_iterator fit= sa.faces_begin(); 
       fit != sa.faces_end(); ++fit){
    sds_.new_face(fit->begin(), fit->end());
  }

  initialize_certificates();
    
  sds_.set_is_building(false);

  for (Sds::Halfedge_iterator it= sds_.halfedges_begin(); 
       it != sds_.halfedges_end(); ++it){
    if (it->curve().is_inside()) check_edge_collapse(it);
  }
  for (Sds::Face_iterator it= sds_.faces_begin(); it != sds_.faces_end(); ++it){
    Halfedge_handle c= it->halfedge();
    do {
      check_edge_face(c);
      c= c->next();
    } while (c != it->halfedge());
  }
  
  audit();
}







void Slice::audit() const {
  sds_.audit();

  /* check for all rules in vertices that if the rule in the vertex
     differs from the rule on the incident edge, they have the same coordinate.
  */
  for (Sds::Vertex_const_iterator vit= sds_.vertices_begin();
       vit != sds_.vertices_end(); ++vit) {
    if (vit->point().is_special()) continue;
    if (!vit->point().is_finite()) {
      CGAL_assertion(sds_.degree(vit) ==2 || sds_.degree(vit) ==3);
    } else if (vit->point().is_sphere_sphere()) {
      CGAL_assertion(sds_.degree(vit) ==4);
    } else {
      CGAL_assertion(sds_.degree(vit) ==3);
    }
    if (vit->point().is_rule_rule()) {
      T::Key x, y;
      Sds::Curve cx, cy;
      Halfedge_const_handle h= vit->halfedge();
      do {
	if (h->curve().is_vertical()) {
	  if (x != T::Key()) {
	    CGAL_assertion(t_.compare_sphere_centers_c(h->curve().key(),
						       x,
						       plane_coordinate(0))
			   == CGAL::EQUAL);
	  } else {
	    cx= h->curve();
	    x= h->curve().key();
	  }
	} else {
	  if (y != T::Key()) {
	    CGAL_assertion(t_.compare_sphere_centers_c(h->curve().key(),
						       y,
						       plane_coordinate(1))
			   == CGAL::EQUAL);
	  } else {
	    cy= h->curve();
	    y= h->curve().key();
	  }
	}
	h= h->opposite()->prev();
      } while (h != vit->halfedge());
      CGAL_assertion(x != T::Key() && y != T::Key());
      CGAL_assertion(t_.compare_sphere_centers_c(vit->point().rule_key(plane_coordinate(0)),
						 x,
						 plane_coordinate(0))
			   == CGAL::EQUAL);
      CGAL_assertion(t_.compare_sphere_centers_c(vit->point().rule_key(plane_coordinate(1)),
						 y,
						 plane_coordinate(1))
			   == CGAL::EQUAL);
    } else if (vit->point().is_sphere_rule()) {
      Sds::Curve k;
      bool hk=false;
      Sds::Curve ko;
      bool hko=false;
      //int dir;
      Sds::Curve c0, c1;
      bool hc0=false, hc1=false;
      
      Halfedge_const_handle h= vit->halfedge();
      do {
	if (h->curve().is_arc()) {
	  if (hc0) {
	    CGAL_assertion(!hc1);
	    c1= h->curve();
	    hc1=true;
	  } else {
	    c0= h->curve();
	    hc0=true;
	  }
	} else {
	  if (!hk) {
	    hk=true;
	    k= h->curve();
	  } else {
	    CGAL_assertion(!hko);
	    hko=true;
	    ko=h->curve();
	  }
	}
	CGAL_assertion(h->vertex() == h->opposite()->prev()->vertex());
	h= h->opposite()->prev();
      } while (h != vit->halfedge());

      CGAL_assertion(c0.key() == c1.key());
      if (c0.other_side() == c1) {

      } else if (c0.key() != k.key()){
	int ed= Sds::Curve::rule_direction(c0, c1);
	Sds::Curve cr= Sds::Curve::make_rule(c0.key(), ed);
	CGAL_assertion(cr.is_vertical() == k.is_vertical());
	if (cr.is_vertical()) {
	  CGAL_assertion(t_.compare_sphere_centers_c(cr.key(),
						     k.key(),
						     plane_coordinate(0))
			 == CGAL::EQUAL);
	} else {
	  CGAL_assertion(t_.compare_sphere_centers_c(cr.key(),
						     k.key(),
						     plane_coordinate(1))
			 == CGAL::EQUAL);
	}
	if (hko) {
	  if (cr.is_vertical()) {
	    CGAL_assertion(t_.compare_sphere_centers_c(cr.key(),
						       ko.key(),
						       plane_coordinate(0))
			 == CGAL::EQUAL);
	  } else {
	    CGAL_assertion(t_.compare_sphere_centers_c(cr.key(),
						       ko.key(),
						       plane_coordinate(1))
			   == CGAL::EQUAL);
	  }
	}
      }
    }
  }
 
}

