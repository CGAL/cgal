#include <CGAL/Arrangement_of_spheres_3/Slice.h>
#include <CGAL/IO/Qt_examiner_viewer_2.h>

#define DPRINT(x) x




void Slice::draw_rz(Qt_examiner_viewer_2 *qtv, NT z) {
  *qtv << CGAL::RED;
  qtv->set_updating_box(true);
  T::Intersect_with_sweep is=tr_.intersect_with_sweep_object(z);
    
    
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
      T::Circle_2 c= is(sphere(hit->curve().key()));
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
  Slice_arrangement sa(spheres_.begin()+3,spheres_.end(), z, inf());
  sds_.clear();
  sds_.set_is_building(true);

  for (Slice_arrangement::Face_iterator fit= sa.faces_begin(); fit != sa.faces_end(); ++fit){
    sds_.new_face(fit->begin(), fit->end());
  }
    
  sds_.set_is_building(false);
}







