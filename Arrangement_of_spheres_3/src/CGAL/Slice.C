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


Slice::DT::Point_2 Slice::display_point_rz(Sds::Point pt, NT z) const {
  T::Sphere_point_3 sp= sphere_point_rz(pt, z);
  //std::cout << "Exact point is " << sp << std::endl;
  return DT::Point_2(sp.approximate_coordinate(Coordinate_index(0)),
		     sp.approximate_coordinate(Coordinate_index(1)));
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
      if (t_.sphere(hit->curve().key()).center().z() > z) {
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
    if (hit->point().first().is_arc()){
      out << hit->point().first().key();
    } else {
      out << hit->point().first();
    }
    out << ":";
    if (hit->point().second().is_arc()){
      out << hit->point().second().key();
    } else {
      out << hit->point().first();
    }
    //out << hit->point().first() << ":" << hit->point().second();
    
    *qtv << CGAL::GRAY;
    *qtv << out.str().c_str();
  }
}









void Slice::set_rz(NT z) {
  Slice_arrangement sa(t_.spheres_begin(),t_.spheres_end(), z, t_.inf());
  sds_.clear();
  sds_.set_is_building(true);

  for (Slice_arrangement::Face_iterator fit= sa.faces_begin(); fit != sa.faces_end(); ++fit){
    sds_.new_face(fit->begin(), fit->end());
  }
    
  sds_.set_is_building(false);
}







