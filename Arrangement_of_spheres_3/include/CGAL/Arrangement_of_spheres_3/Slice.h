#ifndef SLICE_H
#define SLICE_H
#include <CGAL/basic.h>
//#include <CGAL/Halfedge_DS_default.h>
#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_data_structure.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_arrangement.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arrangement_of_spheres_3/Unordered_pair.h>
#include <CGAL/Arrangement_of_spheres_3/Unordered_triple.h>
#include <boost/array.hpp>

#include <CGAL/Arrangement_of_spheres_3/Simulator.h>
#include <CGAL/Kinetic/Simulator_kds_listener.h>

#include <CGAL/IO/Qt_widget.h>

/* invariants
   
*/
struct Slice {
  typedef Slice This;
  typedef Arrangement_of_spheres_traits_3 T;
  typedef Slice_data_structure Sds;
  typedef T::FT NT;
 
  typedef Sds::Halfedge_handle Halfedge_handle;
  typedef Sds::Vertex_handle Vertex_handle;
  typedef Sds::Face_const_handle Face_const_handle;
  typedef Sds::Halfedge_const_handle Halfedge_const_handle;
  typedef Sds::Vertex_const_handle Vertex_const_handle;
  typedef Sds::Face_handle Face_handle;

  typedef CGAL::Simple_cartesian<double> DT;
  typedef Simulator::Event_key Event_key;

  typedef Unordered_pair<T::Key> Intersection_2;
  typedef Unordered_triple<T::Key> Intersection_3;

  struct Rule_event_rep {
    Rule_event_rep(T::Key rk,
		   Rule_direction ri,
		   T::Key ok): rk_(rk), ri_(ri), ok_(ok){}
    bool operator<(const Rule_event_rep &o) const {
      if (rk_ < o.rk_) return true;
      else if (rk_ > o.rk_) return false;
      else if (ri_ < o.ri_) return true;
      else if (ri_ > o.ri_) return false;
      else return ok_ < o.ok_;
    }
    bool operator==(const Rule_event_rep &o) const {
      return rk_== o.rk_ && ri_ == o.ri_ && ok_ == o.ok_;
    }
    const T::Key& rule_key() const {
      return rk_;
    }
    const T::Key& other_key() const {
      return ok_;
    }
    Coordinate_index rule_coordinate() const {
      return other_plane_coordinate(ri_.constant_coordinate());
    }
    Rule_direction rule_direction() const {
      return ri_;
    }

    T::Key rk_;
    Rule_direction ri_;
    T::Key ok_;
  };

  struct On_edge_exception {
    On_edge_exception(Halfedge_handle h): h_(h){}
    Halfedge_handle halfedge_handle() const {
      return h_;
    }
    Halfedge_handle h_;
  };
  struct On_vertex_exception {
    On_vertex_exception(Vertex_handle h): v_(h){}
    Vertex_handle vertex_handle() const {
      return v_;
    }
    Vertex_handle v_;
  };


  class Gui_listener: public Qt_gui::Listener
  {
  public:
    Gui_listener( Qt_gui::Handle &h, Slice *t): Qt_gui::Listener(h), t_(t){}
    virtual void new_notification( Qt_gui::Listener::Notification_type nt) {
      if (nt == Qt_gui::Listener::PICTURE_IS_VALID) {
	t_->draw_rz(Qt_gui::Listener::widget(), Qt_gui::Listener::notifier()->current_time());
      }
    }
  protected:
    Slice *t_;
  };

  friend class Gui_listener;

  typedef ::Simulator Simulator;
  typedef CGAL::Kinetic::Simulator_kds_listener< Simulator::Listener, This> Simulator_listener;
  friend  class CGAL::Kinetic::Simulator_kds_listener< Simulator::Listener, This>;
  


  Slice(T tr);

  Simulator::Handle simulator_handle() {
    return sim_;
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
    KDS functions -----------------------------------------------------
  */
  class Event_base;
  class One_sphere_event;
  class Two_sphere_event;
  class Three_sphere_event;
  class Rule_event;
  class Edge_event;
  
  void process_one_sphere_event(T::Key k);
  void process_intersect_event(T::Key k, T::Key l);
  void process_unintersect_event(T::Key k, T::Key l);
  void process_two_sphere_event(T::Key k, T::Key l, bool first);
  void process_three_sphere_event(T::Key k, T::Key l, T::Key m, bool first);
  void process_rule_collapse_event(T::Key k, Rule_direction r, T::Key o);
  void process_edge_collapse_event(Halfedge_handle h);
  // check for a degeneracy and return true if this handles it.
  void process_degeneracy();
  void rebuild_degenerate(const T::Sphere_point_3 &sp3, Vertex_handle a_vertex);
  bool has_degeneracy() const;
  void initialize_certificates();
  CGAL::Comparison_result compare_concurrent(Event_key k0, Event_key k1) const ;

  
  void audit_one_sphere_event(T::Key k, Event_key ek);
  void audit_two_sphere_event(T::Key k, T::Key l, bool first, Event_key ek);
  void audit_three_sphere_event(T::Key k, T::Key l, T::Key m, bool first, 
				Event_key ek);
  void audit_rule_collapse_event(T::Key k, Rule_direction r, T::Key o, 
				 Event_key ek);
  void audit_edge_collapse_event(Halfedge_handle h, Event_key ek);


  void process_equal_centers_two_sphere_event(T::Key k, T::Key l, bool f);
  void process_aligned_centers_two_sphere_event(T::Key k, T::Key l, bool f,
						CGAL::Comparison_result c[2]);

 // Functions to update certificates ----------------------------------

  void clean_edge(Halfedge_handle h);
  void check_edge_face(Halfedge_handle h);
  void check_edge_pair(Halfedge_handle h0, Halfedge_handle h1);
  void check_merged_faces(Face_handle f, Face_handle g);
  void check_edge_collapse(Halfedge_handle h);
  void check_reduced_face(Face_handle f);

  void set_is_editing(bool tf);



  /*
    Modifiers ----------------------------------------------------------
  */


  /*Face_handle intersect(const T::Sphere_point_3 &fp,
			T::Key ka, T::Key kb);
  Face_handle unintersect(const T::Sphere_point_3 &fp,
  T::Key ka, T::Key kb);*/

  // return the face inside the sphere
  Face_handle insert_sphere(const T::Sphere_point_3 &fp, T::Key k);

  // cur should point out
  T::Key roll_back_rule(const T::Sphere_point_3 &ep,
			Halfedge_handle cur);

  Face_handle erase_sphere(const T::Sphere_point_3 &ep,
			   T::Key k);


  /* Face_handle insert_sphere_on_arc_prep(const T::Sphere_point_3 &ep, 
					T::Key k,
					Halfedge_handle h,
					Halfedge_handle rvs[]);*/
  
  Face_handle insert_sphere_on_rr_prep(const T::Sphere_point_3 &ep, 
				       T::Key k,
				       Vertex_handle v,
				       Vertex_handle rvs[]);
  
  Face_handle insert_sphere_on_rule_prep(const T::Sphere_point_3 &ep, 
					 T::Key k,
					 Halfedge_handle h,
					 Vertex_handle rvs[]);

  /*Face_handle insert_sphere_on_ss_prep(const T::Sphere_point_3 &ep, 
				       T::Key k,
				       Vertex_handle v,
				       Halfedge_handle rvs[]);*/

  Halfedge_handle check_remove_redundant(Halfedge_handle v);

  Face_handle intersect_spheres(const T::Event_point_3 &t,
				Halfedge_handle ha, Halfedge_handle hb);

  Face_handle unintersect_spheres(const T::Event_point_3 &t,
				  Halfedge_handle ha,
				  Halfedge_handle hb);
  
  Face_handle intersect_3_spheres(const T::Event_point_3 &ep,
				  Face_handle f);

  Face_handle collapse_edge(const T::Event_point_3 &ep,
			    Halfedge_handle rule,
			    Halfedge_handle c,
			    Halfedge_handle base);

  Halfedge_handle uncollapse_edge(const T::Event_point_3 &ep,
				  Halfedge_handle rule,
				  Halfedge_handle c,
				  Halfedge_handle base);

  /*
    Replace rule with a rule rotate 90 degres around rule->vertex().
  */
  Halfedge_handle rotate_rule(const T::Event_point_3 &ep,
			      Halfedge_handle rule);

  T::Key debug_new_sphere(T::Sphere_3 s);
  /*
    Display functions----------------------------------------------------
  */


  DT::Point_2 display_point_rz(Sds::Point pt, NT z) const;

  void draw_rz(Qt_examiner_viewer_2 *qtv, NT z);

  void draw_rz(CGAL::Qt_widget *qtv, NT z);

  void draw_events_rz(CGAL::Qt_widget *qtv, NT z);

  void draw_marked_rz(Qt_examiner_viewer_2 *qtv, NT z);

  std::ostream &write(Vertex_const_handle v, std::ostream &out) const;

  std::ostream &write(Halfedge_const_handle v, std::ostream &out) const ; 

  std::ostream &write(Face_const_handle v, std::ostream &out) const;

  /* 
     rational z functions-------------------------------------------------
  */

  void initialize_at(NT z, bool gen_certs=true) ;

  T::Point_2 center_point_rz(T::Key a, T::Key b, NT z) const ;

  T::Circle_2 circle_rz(T::Key a, NT z) const;

  bool intersects_rz(T::Key a, NT z) const;

  T::Sphere_point_3 sphere_point_rz(Sds::Point pt, NT z) const ;

  /*T::Line_3 in_line_rz(Sds::Curve r, NT z) const;
  
  T::Line_3 out_line_rz(Sds::Curve r, NT z) const;*/

  T::Line_3 positive_line_rz(T::Key, Coordinate_index i, NT z) const;
  
  T::Line_3 negative_line_rz(T::Key r, Coordinate_index i, NT z) const;

  /* 
     constructions----------------------------------------------------
  */

  T::Point_2 compute_rule_rule_intersection(T::Key rx, T::Key ry) const ;

  /*
    Search functions ------------------------------------------------
  */


  Face_handle locate_point(const T::Sphere_point_3 &ep);
  
  template <class It>
  Face_handle locate_point(It b, It e, const T::Sphere_point_3 &ep);

  Halfedge_handle shoot_rule(const T::Sphere_point_3& t,
			     Face_handle f,
			     const T::Sphere_point_3& source,
			     Rule_direction ruledir) ;

  Halfedge_handle 
  find_rule_vertex(const T::Sphere_point_3 &t, 
		   Face_handle f,
		   Sds::Curve rule);
  
  /* 
     predictes--------------------------------------------------------
  */


  //enum Location {L_BIT=1, R_BIT=2, T_BIT=4, B_BIT=8, IN_BIT=16,  };

  int sphere_location(const T::Sphere_point_3 &sp,
		      T::Key sphere) const ;
  
  bool behind_arc(const T::Sphere_point_3 &ep,
		  Sds::Curve arc, int location) const;

  void point_sphere_orientation(const T::Sphere_point_3 &time,
				T::Key sphere,
				std::vector<int> &locations) const;

  bool locate_point_check_face(const T::Sphere_point_3 &z,
			       Face_const_handle it,
			       std::vector<int> &locations) const ;

  bool locate_point_check_face_arcs(const T::Sphere_point_3 &z,
				    Face_const_handle f,
				    std::vector<int> &locations) const ;

  bool locate_point_check_face_vertices(const T::Sphere_point_3 &ep,
					Face_const_handle it) const;

  CGAL::Comparison_result debug_rule_shoot_answer(const T::Sphere_point_3 &t,
						  const T::Sphere_point_3 &source,
						  Rule_direction ruledir,
						  Sds::Point pt,
						  bool &exact) const ;

  
  
  static void debug_rule_shoot_check(CGAL::Comparison_result check, 
				     CGAL::Comparison_result computed,
				     bool exact) ;



  CGAL::Comparison_result rule_shoot_compare_SR(const T::Sphere_point_3 &t,
						const T::Sphere_point_3 &pt,
						Rule_direction ruledir,
						Sds::Curve arc,
						T::Key orule,
						Sds::Point debug_pt,
						bool arc_above) const ;



  CGAL::Comparison_result rule_shoot_compare_SS(const T::Sphere_point_3 &t,
						const T::Sphere_point_3 &pt,
						Rule_direction ruledir,
						Sds::Point pt) const ;
 

  // return comparison of separator to intersection point on the C coordinate
  // i.e. SMALLER if the separator is SMALLER than the intersection point
  CGAL::Comparison_result rule_shoot_edge_vertex(const T::Sphere_point_3 &t,
						 const T::Sphere_point_3 &pt,
						 Rule_direction rule,
						 Sds::Curve hp,
						 Sds::Point p,
						 Sds::Curve hn) const ;
  
  bool rule_shoot_compare_if_rational(const T::Sphere_point_3 &pt, 
				      const T::Sphere_point_3 &t, 
				      Rule_direction rule,
				      Sds::Point pt,
				      CGAL::Comparison_result &ret) const ;
  /*bool rule_shoot_compare_if_rational_arc(const T::Sphere_point_3 &z,
					  Sds::Curve rule,
					  Sds::Curve a,
					  CGAL::Comparison_result &ret) const;*/
  

 


  // Debug Functions ---------------------------------------------------

  //Face_handle locate_point(const T::Sphere_point_3 &ep);

  void audit() const;
  void audit_events() const;

  T &traits() {
    return t_;
  }


  void set_gui(Qt_gui::Handle qt);

private:
  typedef std::pair<Event_key,Event_key> Event_key_pair;

  struct Handle_compare{
    template <class H>
    bool operator()(H a, H b) const {
      return &*a < &*b;
    }
  };

  /*struct Rule_events {
    Rule_events(){}

    Simulator::Event_key collapse() const {
      return collapse_;
    }
    Simulator::Event_key slide(bool pos) {
      if (pos) return slide_[1];
      else return slide_[0];
    }
    void set_collapse(Simulator::Event_key k) {
      collapse_=k;
    }
    void set_slide(bool pos, Simulator::Event_key k) {
      if (pos) slide_[1]=k;
      else slide_[0]=k;
    }
    
    Simulator::Event_key collapse_, slide_[2];
    };*/

  T t_;
  Slice_data_structure sds_;
  

  typedef std::map<Intersection_2, Event_key_pair > Intersections_2;
  Intersections_2 intersections_2_;
  typedef std::map<Intersection_3, Event_key_pair > Intersections_3;
  Intersections_3 intersections_3_;
  typedef std::map<Rule_event_rep, Event_key_pair > Rule_events;
  Rule_events rule_events_;
  /*std::map<Unordered_pair<T::Key>, Face_handle> pair_map_;
    std::map<Unordered_triple<T::Key>, Face_handle> triple_map_;*/

  std::set<Face_const_handle, Handle_compare> marked_faces_;
  std::set<Vertex_const_handle, Handle_compare> marked_vertices_;
  std::set<Halfedge_const_handle, Handle_compare> marked_edges_;

  Simulator::Handle sim_;
  Simulator_listener siml_;
  Gui_listener *guil_;
};

#endif
