#ifndef CGAL_AOS3_DSITEMSEVENT_VISITOR_H
#define CGAL_AOS3_DSITEMSEVENT_VISITOR_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Rule_direction.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_vertex.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve.h>


CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
struct Cross_section_halfedgeDS_items_2 {
  CGAL_AOS3_TRAITS;
public:
  typedef Combinatorial_vertex Point;
  typedef Combinatorial_curve Curve;
  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Rule_direction Rule_direction;
 
  typedef CGAL_AOS3_TYPENAME Traits::Event_key Event_key;

  template < class Refs, class Traits>
  struct Vertex_wrapper {
     
    struct Vertex: public CGAL::HalfedgeDS_vertex_base< Refs, CGAL::Tag_true, Point>{
    };
  };
  template < class Refs, class Traits>
  struct Halfedge_wrapper {
    struct Halfedge: public CGAL::HalfedgeDS_halfedge_base< Refs> {
      Halfedge(){}
      Curve curve() const {
	return pt_;
      }
      Curve& curve() {
	return pt_;
      }
      void set_curve(Curve pt) {
	pt_=pt;
      }
      Event_key event() const {
	return ev_;
      }
      void set_event(Event_key ev) {
	ev_=ev;
      }

      Event_key ev_;
      Curve pt_;
    };
  };
  template < class Refs, class Traits>
  struct Face_wrapper {
    // maybe want pointers to chains, maybe not
    typedef CGAL::HalfedgeDS_face_base< Refs> Face;
  };
};
  
CGAL_AOS3_END_INTERNAL_NAMESPACE
#endif
