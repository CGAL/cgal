#ifndef CGAL_AOS3_EVENT_VISITOR_H
#define CGAL_AOS3_EVENT_VISITOR_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_halfedgeDS_items_2.h>
//#include <CGAL/Arrangement_of_spheres_3/Event_processor.h>
#include <set>
#include <CGAL/box_intersection_d.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Event_processor ;

CGAL_AOS3_TEMPLATE
class Combinatorial_cross_section ;

/*
  Event_visitor gets types from CCS-- why bother?
  Event_handler manipulates CCS
  CCS calls Event_visitor
 */

CGAL_AOS3_TEMPLATE 
struct Event_visitor: boost::noncopyable {
  CGAL_AOS3_TRAITS;
  typedef Cross_section_halfedgeDS_items_2 CGAL_AOS3_TARG HDS_items;
  typedef Event_processor CGAL_AOS3_TARG EP;
public:
  typedef CGAL::HalfedgeDS_default<int, HDS_items> HDS;

  typedef CGAL_AOS3_TYPENAME HDS::Halfedge_handle Halfedge_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Halfedge_const_handle Halfedge_const_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Vertex_handle Vertex_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Vertex_const_handle Vertex_const_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Face_handle Face_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Face_const_handle Face_const_handle;
  typedef CGAL_AOS3_TYPENAME Traits::Event_key Event_key;
  typedef CGAL_AOS3_TYPENAME HDS_items::Curve Curve;
  typedef CGAL_AOS3_TYPENAME HDS_items::Point Point;
  typedef CGAL_AOS3_TYPENAME Curve::Key Sphere_3_key;
  typedef CGAL_AOS3_TYPENAME Traits::Event_pair Event_pair;
 
  typedef CGAL_AOS3_TYPENAME Traits::Simulator Simulator;

  Event_visitor():has_traits_(false),  j_(NULL){}

  void set_traits(Traits &tr);
  void set_cross_section(Combinatorial_cross_section CGAL_AOS3_TARG &cs);

  
  
  void on_new_edge(Halfedge_handle h);
  void on_delete_edge(Halfedge_handle h);
  void on_change_edge(Halfedge_handle h);
  void on_merge_faces(Halfedge_handle h);

  void audit(Halfedge_const_handle h) const;
  void set_start_time(CGAL_AOS3_TYPENAME Traits::FT st) {
    start_time_=st;
    has_start_time_=true;
  }
  
  CGAL_AOS3_TYPENAME Simulator::Handle simulator() {
    return sim_;
  }

  CGAL_AOS3_TYPENAME Simulator::Const_handle simulator() const {
    return sim_;
  }

  CGAL_CONST_ITERATOR(Free_event,free_event,
		      CGAL_AOS3_TYPENAME std::set<Event_key>::const_iterator,
		      return free_events_.begin(),
		      return free_events_.end());
  
  void delete_free_event() {
    free_events_.erase(sim_->current_event());
  }

private:


typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,std::pair<Sphere_3_key, Sphere_3_key>,
						      Box_intersection_d::ID_EXPLICIT> Box;

  struct Report_pairs {
    Report_pairs(Event_visitor CGAL_AOS3_TARG *v): v_(v){}
    template <class B>
    void operator()(const B& a, const B& b) {
      v_->process_pair(a,b);
    }
    Event_visitor CGAL_AOS3_TARG *v_;
  };

  struct Report_triples {
    Report_triples(Event_visitor CGAL_AOS3_TARG *v): v_(v){}
    template <class B>
    void operator()(const B& a, const B& b) {
      CGAL_precondition(b.handle().first != b.handle().second);
      CGAL_precondition(a.handle().first == a.handle().second);
      if (b.handle().first != a.handle().first 
	  && b.handle().first != a.handle().second
	  && a.handle().first < b.handle().first 
	  && a.handle().first < b.handle().second) {
	v_->process_triple(a.handle().first,
			   b.handle().first,
			   b.handle().second);
      }
    }
    Event_visitor CGAL_AOS3_TARG *v_;
  };

  bool should_have_certificate(Halfedge_const_handle h) const ;

  void initialize();

  void new_event(Halfedge_handle h);
    
  void handle_edge_face(Halfedge_handle h, Face_handle f) ;


  std::ostream &write(Vertex_const_handle h, 
					  std::ostream &out) const {
  out << h->point();
  return out;
}


  
  std::ostream &write(Halfedge_const_handle h, 
		      std::ostream &out) const {
  if (h == Halfedge_const_handle()) out << "NULL";
  else {
    out << h->opposite()->vertex()->point() << " -- " << h->curve()
	<< " -- " << h->vertex()->point();
  }
  return out;
}


  void set_event(Halfedge_handle h, Event_key k);
  void process_pair(const Box &a, const Box &b);
  void process_pair(Sphere_3_key a, 
		    Sphere_3_key b);
  
  void process_triple(Sphere_3_key a, 
		      Sphere_3_key b,
		      Sphere_3_key c);
    
  //void process_triple(Halfedge_handle h);
  
  
  //std::set< Sphere_key_upair> pairs_;
  // std::set< Sphere_key_utriple> triples_;
  Traits tr_;
  bool has_traits_;
  CGAL_AOS3_TYPENAME Traits::FT start_time_;
  bool has_start_time_;
  std::set<Event_key> free_events_;

  EP * j_;
  CGAL_AOS3_TYPENAME Simulator::Handle sim_;


  
  std::vector<Box> intersections_;
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Event_visitor_impl.h>
#endif

#endif
