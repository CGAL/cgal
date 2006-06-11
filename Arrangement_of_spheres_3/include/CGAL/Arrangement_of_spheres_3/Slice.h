#ifndef SLICE_H
#define SLICE_H
#include <CGAL/basic.h>
//#include <CGAL/Halfedge_DS_default.h>
#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_data_structure.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_arrangement.h>
#include <CGAL/Arrangement_of_spheres_3/utilities.h>

/* invariants
   
*/
struct Slice {
  typedef Arrangement_of_spheres_traits_3 T;
  typedef Slice_data_structure Sds;
  typedef T::NT NT;
  typedef Sds::Face_handle Face_handle;
  typedef Sds::Halfedge_const_handle Halfedge_const_handle;
  typedef Sds::Face_const_handle Face_const_handle;
  typedef Sds::Vertex_const_handle Vertex_const_handle;
  typedef Sds::Face_const_iterator Face_const_iterator;

  struct On_edge_exception {
    On_edge_exception(Halfedge_const_handle h): h_(h){}
    Halfedge_const_handle halfedge_handle() const {
      return h_;
    }
    Halfedge_const_handle h_;
  };
  struct On_vertex_exception {
    On_vertex_exception(Vertex_const_handle h): v_(h){}
    Vertex_const_handle vertex_handle() const {
      return v_;
    }
    Vertex_const_handle v_;
  };


  template <class It> 
  Slice(It bs, It es, NT inf): spheres_(bs, es), inf_(inf){
  }
 

  
 


 

  void new_marked_face(Face_const_handle f) {
    marked_faces_.insert(f);
  }

  void marked_faces_clear() {
    marked_faces_.clear();
  }
  void new_marked_vertex(Vertex_const_handle f) {
    marked_vertices_.insert(f);
  }

  void marked_vertices_clear() {
    marked_vertices_.clear();
  }
  void new_marked_edge(Halfedge_const_handle f) {
    marked_edges_.insert(f);
  }

  void marked_edges_clear() {
    marked_edges_.clear();
  }


  std::ostream &write(Vertex_const_handle v, std::ostream &out) const {
    out << v->point();
    return out;
  }
  std::ostream &write(Halfedge_const_handle v, std::ostream &out) const {
    out << v->curve();
    return out;
  }
  std::ostream &write(Face_const_handle v, std::ostream &out) const {
    Halfedge_const_handle h=v->halfedge();
    do {
      out << h->curve() << "--" << h->vertex()->point() << "--";
      h=h->next();
    } while (h != v->halfedge());
    return out;
  }

  //Face_const_handle locate_point(int ind) const;

  Face_const_handle locate_point(int ind, T::Event_point_3 ep) const;

  Face_const_handle locate_point(T::Event_point_3 ep) const;


  Halfedge_const_handle shoot_rule(T::Event_point_3 source, Face_const_handle f,
				   Sds::Curve rule) const ;
  
  Halfedge_const_handle shoot_rule(T::Event_point_3 source, Face_const_handle f,
				   int type) const ;


  bool build_intersection_events(int a, int b, T::Event_point_3 &begin, T::Event_point_3 &end) const {
    // build equipower plane
    // intersect it with forward plane (check which is needed)
    // build events
    
    // check if events are valid
    CGAL_assertion(0);
    return false;
  }

  bool intersection_event(Halfedge_const_handle &a, Halfedge_const_handle &b, 
			  T::Event_point_3 &begin, T::Event_point_3 &end) const {
    CGAL_precondition(!a->curve().is_rule());
    CGAL_precondition(!b->curve().is_rule());
    if (a->curve().is_inside() && b->curve().is_inside()) {
      return false;
    }
    if (b->curve().is_inside()) {
      std::swap(a,b);
    }
    // b is not inside
    if (a->curve().is_inside()) {
      if (a->curve().quadrant() == b->curve().quadrant()){
	return build_intersection_events(a->curve().index(), b->curve().index(), 
					 begin, end);
      } else return false;
    } else {
      if (a->curve().quadrant() & b->curve().quadrant()) {
	return false;
      } else return build_intersection_events(a->curve().index(),
					      b->curve().index(), begin, end);
    }
  }

  // in a degeneracy this could result in a circle being formed
  bool rule_collapse_event(Halfedge_const_handle a, T::Event_point_3 &begin, 
			   T::Event_point_3 &end) const {
    return false;
  }

  bool face_collapse_event(Face_const_handle f, T::Event_point_3 &t) const {
    // make sure face has three arc edges (two arc edges or just one circle will be handled separately)
    // make two equipower planes
    // intersect them
    // they must intersect
    // check line against one sphere
    return false;
  }


 
 

 

 

 
 
  



  // rational z functions-------------------------------------------------

  void draw_rz(Qt_examiner_viewer_2 *qtv, NT z) ;

  void set_rz(NT z) ;

  T::Point_2 center_point_rz(int a, int b, NT z) const ;

  T::Point_2 display_point_rz(Sds::Point pt, NT z) const;

  T::Sphere_point_3 sphere_point_rz(Sds::Point pt, NT z) const ;

  /* 
     constructions----------------------------------------------------
   */

 // the point described by the vertex (a,b) should be on the positive side
  T::Plane_3 separating_plane(int a, int b) const ;

  // point from the second to the first
  T::Plane_3 equipower_plane(int a, int b) const ;

  NT disc(int i) const;

  NT center_c(int ind, int C) const;

  T::Point_3 center(int ind) const;

  T::Point_2 compute_rule_rule_intersection(Sds::Curve ra, Sds::Curve rb) const ;

  T::Sphere_3 sphere(int ind) const;


  T::Line_3 in_line(Sds::Curve r, NT z) const;
  
  T::Line_3 out_line(Sds::Curve r, NT z) const;

  T::Event_point_3 sphere_start(int s) const;

 
  /*
    Objects
  */

  // extremal points-- lines point forward or backwards
  // general sphere points

  /* 
     predictes--------------------------------------------------------
   */

  bool intersects_rule(int sphere, int rule_sphere, int C) const;
  
  //CGAL::Comparison_result compare_sphere_to_plane(int s, NT c, int C) const;

  CGAL::Comparison_result compare_equipower_point_to_rule(int sphere0,
							  int sphere1,
							  int rule_sphere,
							  int C) const;

  CGAL::Comparison_result compare_sphere_center_to_rule(int sphere, 
							int rule_sphere,
							int C) const;

  

  /*CGAL::Comparison_result compare_equipower_point_to_plane(int a, int b, 
    NT coord, int C) const;*/

  CGAL::Sign sign_of_separating_plane_normal_c(int sphere0, int sphere1,
					     int C) const ;


  CGAL::Sign sign_of_equipower_plane_normal_c(int sphere0, int spher1, 
					      int C) const;


  CGAL::Oriented_side oriented_side_of_equipower_plane(int sphere_0, int sphere_1,
						       const T::Sphere_point_3 &s) const;

  CGAL::Oriented_side oriented_side_of_center_plane(int sphere_0, int sphere_1, 
						    int sphere_center) const ;

 
  CGAL::Comparison_result compare_sphere_centers_c(int a, int b, int C) const;







  /*
    High level predicates
  */

  //enum Location {L_BIT=1, R_BIT=2, T_BIT=4, B_BIT=8, IN_BIT=16,  };

  int sphere_location(int locate_point, int sphere) const ;
  
  bool behind_arc(T::Event_point_3 ep, int ind, Sds::Curve arc, int location) const;


  void point_sphere_orientation(int point,
				int sphere,
				std::vector<int> &locations) const;

  bool locate_point_check_face(Face_const_handle it, 
			       int ind,
			       std::vector<int> &locations) const ;

  bool locate_point_check_face_arcs(T::Event_point_3 ep,
				    int ind,
				    Face_const_handle f,
				    std::vector<int> &locations) const ;

  bool locate_point_check_face_vertices(T::Event_point_3 ep,
					Face_const_handle it) const;

  CGAL::Comparison_result debug_rule_shoot_answer(T::Event_point_3 ep, 
						  Sds::Curve rule,
						  Sds::Curve p,
						  Sds::Curve n,
						  bool &exact) const ;

  
  
  static void debug_rule_shoot_check(CGAL::Comparison_result check, 
				     CGAL::Comparison_result computed,
				     bool exact) ;



  CGAL::Comparison_result rule_shoot_compare_SR(T::Event_point_3 ep, 
						Sds::Curve srule,
						Sds::Curve arc,
						Sds::Curve orule,
						bool arc_above) const ;



 CGAL::Comparison_result rule_shoot_compare_SS(T::Event_point_3 ep, 
					       Sds::Curve srule,
					       Sds::Curve arc0,
					       Sds::Point pt,
					       Sds::Curve arc1) const ;
 

  // return comparison of separator to intersection point on the C coordinate
  // i.e. SMALLER if the separator is SMALLER than the intersection point
  CGAL::Comparison_result rule_shoot_edge_vertex(T::Event_point_3 ep, 
						 Sds::Curve rule,
						 Sds::Curve hp,
						 Sds::Point p,
						 Sds::Curve hn) const ;

  bool rule_shoot_compare_if_rational(T::Event_point_3 ep, 
				      Sds::Curve rule,
				      Sds::Curve a,
				      Sds::Curve b,
				      CGAL::Comparison_result &ret) const ;
  bool rule_shoot_compare_if_rational_arc(T::Event_point_3 ep,
					  Sds::Curve rule,
					  Sds::Curve a,
					  CGAL::Comparison_result &ret) const;
  
  // compare along C where the intersection point of the two rules
  // (defined by Coord C of ep) is relative to the sphere of the arc
  CGAL::Bounded_side rule_shoot_source_side(T::Event_point_3 ep,
					    Sds::Curve orule,
					    int arc_index,
					    int C) const;

private:

  struct Handle_compare{
    template <class H>
    bool operator()(H a, H b) const {
      return &*a < &*b;
    }
  };


  struct Temp_point {
    Temp_point(std::vector<T::Sphere_3>  &ss, T::Sphere_3 s): ss_(ss){
      ss_.push_back(s);
    }
    ~Temp_point(){
      ss_.pop_back();
    }
    int index() const {
      return ss_.size()-1;
    }
    std::vector<T::Sphere_3> &ss_;
  };

  //  HDS hds_;
  // ick, this is to handle location of points which are not already there
  // -1, -2 are bl, tr inf corners
  mutable std::vector<T::Sphere_3> spheres_;
  std::vector<int> errors_;
  NT inf_;
  Slice_data_structure sds_;
  T tr_;
  std::set<Face_const_handle, Handle_compare> marked_faces_;
  std::set<Vertex_const_handle, Handle_compare> marked_vertices_;
  std::set<Halfedge_const_handle, Handle_compare> marked_edges_;
};

#endif
