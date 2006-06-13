#ifndef SLICE_H
#define SLICE_H
#include <CGAL/basic.h>
//#include <CGAL/Halfedge_DS_default.h>
#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_data_structure.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_arrangement.h>
#include <CGAL/Arrangement_of_spheres_3/utilities.h>
#include <CGAL/Simple_cartesian.h>

/* invariants
   
*/
struct Slice {
  typedef Arrangement_of_spheres_traits_3 T;
  typedef Slice_data_structure Sds;
  typedef T::FT NT;
  typedef Sds::Face_handle Face_handle;
  typedef Sds::Halfedge_const_handle Halfedge_const_handle;
  typedef Sds::Face_const_handle Face_const_handle;
  typedef Sds::Vertex_const_handle Vertex_const_handle;
  typedef Sds::Face_const_iterator Face_const_iterator;

  typedef CGAL::Simple_cartesian<double> DT;

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
  Slice(It bs, It es): t_(bs, es){
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


 

  //Face_const_handle locate_point(int ind) const;

 

  /*bool intersection_event(Halfedge_const_handle &a, Halfedge_const_handle &b, 
			  T::Sphere_point_3 &begin, T::Sphere_point_3 &end) const {
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
	return build_intersection_events(a->curve().key(), b->curve().key(), 
					 begin, end);
      } else return false;
    } else {
      if (a->curve().quadrant() & b->curve().quadrant()) {
	return false;
      } else return build_intersection_events(a->curve().key(),
					      b->curve().key(), begin, end);
    }
    }*/

  // in a degeneracy this could result in a circle being formed
  /*bool rule_collapse_event(Halfedge_const_handle a, T::Sphere_point_3 &begin, 
			   T::Sphere_point_3 &end) const {
    return false;
    }*/

  /* bool face_collapse_event(Face_const_handle f, T::Sphere_point_3 &t) const {
    // make sure face has three arc edges (two arc edges or just one circle will be handled separately)
    // make two equipower planes
    // intersect them
    // they must intersect
    // check line against one sphere
    return false;
    }*/



  /*
    Modifiers ----------------------------------------------------------
  */

  // return the face inside the sphere
  Face_const_handle insert_sphere(const T::Sphere_point_3 &fp, T::Key k);




  /*
    Display functions----------------------------------------------------
  */


  DT::Point_2 display_point_rz(Sds::Point pt, NT z) const;

  void draw_rz(Qt_examiner_viewer_2 *qtv, NT z) ;

  std::ostream &write(Vertex_const_handle v, std::ostream &out) const;

  std::ostream &write(Halfedge_const_handle v, std::ostream &out) const ; 

  std::ostream &write(Face_const_handle v, std::ostream &out) const;

  /* 
     rational z functions-------------------------------------------------
  */

  void set_rz(NT z) ;

  T::Point_2 center_point_rz(T::Key a, T::Key b, NT z) const ;

  T::Circle_2 circle_rz(T::Key a, NT z) const;

  bool intersects_rz(T::Key a, NT z) const;

  T::Sphere_point_3 sphere_point_rz(Sds::Point pt, NT z) const ;

  /* 
     constructions----------------------------------------------------
   */

  T::Point_2 compute_rule_rule_intersection(Sds::Curve ra, Sds::Curve rb) const ;

  T::Line_3 in_line(Sds::Curve r, NT z) const;
  
  T::Line_3 out_line(Sds::Curve r, NT z) const;


  /*
    Search functions ------------------------------------------------
  */


  Face_const_handle locate_point(const T::Sphere_point_3 &ep, T::Key ind) const;

  Halfedge_const_handle shoot_rule(const T::Sphere_point_3& source,
				   Face_const_handle f,
				   Sds::Curve rule) const ;

  /* 
     predictes--------------------------------------------------------
   */


  //enum Location {L_BIT=1, R_BIT=2, T_BIT=4, B_BIT=8, IN_BIT=16,  };

  int sphere_location(const T::Sphere_point_3 &sp,
		      T::Key locate_point,
		      T::Key sphere) const ;
  
  bool behind_arc(const T::Sphere_point_3 &ep, T::Key ind,
		  Sds::Curve arc, int location) const;


  void point_sphere_orientation(const T::Sphere_point_3 &time,
				T::Key point,
				T::Key sphere,
				std::vector<int> &locations) const;

  bool locate_point_check_face(const T::Sphere_point_3 &z,
			       Face_const_handle it, 
			       T::Key ind,
			       std::vector<int> &locations) const ;

  bool locate_point_check_face_arcs(const T::Sphere_point_3 &z,
				    T::Key ind,
				    Face_const_handle f,
				    std::vector<int> &locations) const ;

  bool locate_point_check_face_vertices(const T::Sphere_point_3 &ep,
					T::Key index,
					Face_const_handle it) const;

  CGAL::Comparison_result debug_rule_shoot_answer(const T::Sphere_point_3 &z,
						  Sds::Curve rule,
						  Sds::Curve p,
						  Sds::Curve n,
						  bool &exact) const ;

  
  
  static void debug_rule_shoot_check(CGAL::Comparison_result check, 
				     CGAL::Comparison_result computed,
				     bool exact) ;



  CGAL::Comparison_result rule_shoot_compare_SR(const T::Sphere_point_3 &z,
						Sds::Curve srule,
						Sds::Curve arc,
						Sds::Curve orule,
						bool arc_above) const ;



 CGAL::Comparison_result rule_shoot_compare_SS(const T::Sphere_point_3 &z,
					       Sds::Curve srule,
					       Sds::Curve arc0,
					       Sds::Point pt,
					       Sds::Curve arc1) const ;
 

  // return comparison of separator to intersection point on the C coordinate
  // i.e. SMALLER if the separator is SMALLER than the intersection point
  CGAL::Comparison_result rule_shoot_edge_vertex(const T::Sphere_point_3 &z,
						 Sds::Curve rule,
						 Sds::Curve hp,
						 Sds::Point p,
						 Sds::Curve hn) const ;

  bool rule_shoot_compare_if_rational(const T::Sphere_point_3 &z, 
				      Sds::Curve rule,
				      Sds::Curve a,
				      Sds::Curve b,
				      CGAL::Comparison_result &ret) const ;
  bool rule_shoot_compare_if_rational_arc(const T::Sphere_point_3 &z,
					  Sds::Curve rule,
					  Sds::Curve a,
					  CGAL::Comparison_result &ret) const;
  


  // Debug Functions ---------------------------------------------------
  Face_const_handle locate_point(const T::Sphere_point_3 &ep) const;

  Halfedge_const_handle shoot_rule(const T::Sphere_point_3 &source, 
				   Face_const_handle f,
				   int type) const ;
  
  Face_const_handle insert_sphere(const T::Sphere_point_3 &ep) const;

private:

  struct Handle_compare{
    template <class H>
    bool operator()(H a, H b) const {
      return &*a < &*b;
    }
  };




  

  /*struct Temp_point {
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
    };*/

  //  HDS hds_;
  // ick, this is to handle location of points which are not already there
  // -1, -2 are bl, tr inf corners
  T t_;
  std::vector<int> errors_;
  Slice_data_structure sds_;
  std::set<Face_const_handle, Handle_compare> marked_faces_;
  std::set<Vertex_const_handle, Handle_compare> marked_vertices_;
  std::set<Halfedge_const_handle, Handle_compare> marked_edges_;
};

#endif
