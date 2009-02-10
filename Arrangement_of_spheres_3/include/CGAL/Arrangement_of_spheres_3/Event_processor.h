#ifndef CGAL_AOS3_JANITOR_H
#define CGAL_AOS3_JANITOR_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>
#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section_rules.h>
#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section_location.h>
#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section_insertion.h>
#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section_removal.h>
#include <CGAL/Arrangement_of_spheres_3/Degenerate_cross_section.h>
#include <CGAL/Kinetic/Event_base.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE


/*
  Moving points are:
  arc/arc AA
  arc/rule AR
  extremum E

  stationary points are:
  rule/rule RR

  events are:
  AA-AA -intersect 3
  AA-AR
  AA-E -- don't deschedule
  AR RR
  AR-AR
  E RR


  Finally:
  AAR
  RAR
  ARR

  AAE

 */
CGAL_AOS3_TEMPLATE
class Event_processor : public CGAL::Kinetic::Ref_counted<Event_processor CGAL_AOS3_TARG> {
  CGAL_AOS3_TRAITS;
  typedef Event_processor CGAL_AOS3_TARG This;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CCS;
  typedef Irrational_cross_section_insertion CGAL_AOS3_TARG ICSI;
  typedef Degenerate_cross_section CGAL_AOS3_TARG DCS;
  typedef Irrational_cross_section_removal CGAL_AOS3_TARG ICSD;
  typedef Irrational_cross_section_location CGAL_AOS3_TARG ICSL;
  typedef Irrational_cross_section_rules CGAL_AOS3_TARG ICSR;
  typedef CGAL_AOS3_TYPENAME Traits::Sphere_3_key Sphere_3_key;
 
  typedef CGAL_AOS3_TYPENAME Traits::Degeneracy_exception Degeneracy;
  typedef CGAL_AOS3_TYPENAME Traits::Event_point_3 Event_point_3;

public:
 
  typedef CGAL_AOS3_TYPENAME CCS::Halfedge_handle Halfedge_handle;
  typedef CGAL_AOS3_TYPENAME CCS::Vertex_handle Vertex_handle;
  typedef CGAL_AOS3_TYPENAME CCS::Face_handle Face_handle;
  typedef CGAL_AOS3_TYPENAME CCS::Point Point;
  typedef CGAL_AOS3_TYPENAME CCS::Curve Curve;
  
 struct Event_base: public CGAL::Kinetic::Event_base<This*> {
    virtual bool is_edge_event() const {return false;}
    virtual std::ostream &write(std::ostream&out) const{return out;};
    virtual ~Event_base(){}
    Event_base(This* t): CGAL::Kinetic::Event_base<This*>(t) {}
  };

#define CGAL_AOS3_EVENT_1K(UCname, lcname)				\
  struct UCname ##_event: Event_base {					\
    UCname ##_event(This* h, Sphere_3_key k): Event_base(h), k_(k){}	\
    void process() {							\
      Event_base::kds()->delete_free_event();				\
      Event_base::kds()->lcname(k_);					\
    }									\
    std::ostream &write(std::ostream &out) const {			\
      out << #UCname << ":"  << k_;					\
      return out;							\
    }									\
    template <class K>							\
    CGAL::Comparison_result compare_concurrent(K a, K b) const {	\
      return Event_base::kds()->compare_concurrent(a,b);		\
    }									\
    Sphere_3_key k_;							\
  };									\
 void lcname(Sphere_3_key k)

#define CGAL_AOS3_EVENT_2K(UCname, lcname)				\
  struct UCname ##_event: Event_base {					\
    UCname ##_event(This* h, Sphere_3_key k, Sphere_3_key l):		\
    Event_base(h), k_(k), l_(l){}					\
    void process() {							\
      Event_base::kds()->delete_free_event();				\
      Event_base::kds()->lcname(k_, l_);				\
    }									\
    std::ostream &write(std::ostream &out) const {			\
      out << #UCname << ":" << k_ << ":" << l_;				\
      return out;							\
    }									\
    template <class K>							\
    CGAL::Comparison_result compare_concurrent(K a, K b) const {	\
      return Event_base::kds()->compare_concurrent(a,b);		\
    }									\
    Sphere_3_key k_, l_;						\
  };									\
 void lcname(Sphere_3_key k, Sphere_3_key l)

#define CGAL_AOS3_EVENT_3K(UCname, lcname)				\
  struct UCname ##_event: Event_base {					\
    UCname ##_event(This* h, Sphere_3_key k, Sphere_3_key l, Sphere_3_key m): \
    Event_base(h), k_(k), l_(l), m_(m){}				\
    void process() {							\
      Event_base::kds()->delete_free_event();				\
      Event_base::kds()->lcname(k_, l_, m_);				\
    }									\
    std::ostream &write(std::ostream &out) const {			\
      out << #UCname << ":" <<  k_ << ":" << l_ << ":" << m_;		\
      return out;							\
    }									\
    template <class K>							\
    CGAL::Comparison_result compare_concurrent(K a, K b) const {	\
      return Event_base::kds()->compare_concurrent(a,b);		\
    }									\
    Sphere_3_key k_, l_, m_;						\
  };									\
 void lcname(Sphere_3_key k, Sphere_3_key l, Sphere_3_key m)

#define CGAL_AOS3_EVENT_1H(UCname, lcname, stype, etype)		\
  struct UCname ##_event: Event_base {					\
    UCname ##_event(This* h, Halfedge_handle hh):			\
    Event_base(h), h_(hh){						\
      CGAL_assertion(hh->opposite()->vertex()->point().stype());	\
      CGAL_assertion(hh->vertex()->point().etype());			\
    }									\
    void process() {							\
      Event_base::kds()->lcname(h_);					\
    }									\
    std::ostream &write(std::ostream &out) const {			\
      out << #UCname << h_->opposite()->vertex()->point()		\
	  << "--" << h_->curve() << "--" << h_->vertex()->point();	\
      return out;							\
    }									\
    template <class K>							\
    CGAL::Comparison_result compare_concurrent(K a, K b) const {	\
      return Event_base::kds()->compare_concurrent(a,b);		\
    }									\
    template <class K>							\
    void audit(K k) const {						\
      CGAL_assertion(h_->event()==k);					\
    }									\
    bool is_edge_event() const {return true;}				\
    Halfedge_handle h_;							\
  };									\
void lcname(Halfedge_handle h)

#define CGAL_AOS3_EVENT_2K1I(UCname, lcname)				\
  struct UCname ##_event: Event_base {					\
    UCname ##_event(This* h, Sphere_3_key k, Sphere_3_key l, int i):	\
    Event_base(h), k_(k), l_(l), i_(i){}				\
    void process() {							\
      Event_base::kds()->delete_free_event();				\
      Event_base::kds()->lcname(k_, l_, i_);				\
    }									\
    std::ostream &write(std::ostream &out) const {			\
      out << #UCname << ":" << k_ << ":" << l_ << "-" << i_;	\
      return out;							\
    }									\
    template <class K>							\
    CGAL::Comparison_result compare_concurrent(K a, K b) const {	\
      return Event_base::kds()->compare_concurrent(a,b);		\
    }									\
    Sphere_3_key k_, l_;						\
    int i_;								\
  };									\
  void lcname(Sphere_3_key k, Sphere_3_key l, int i)

 Comparison_result compare_concurrent(CGAL_AOS3_TYPENAME Traits::Simulator::Event_key a,
				      CGAL_AOS3_TYPENAME Traits::Simulator::Event_key b) const {
   for (int i=0; i< 2; ++i) {
     Comparison_result cr= cs_.visitor().simulator()->event_time(a).compare_c(cs_.visitor().simulator()->event_time(b), plane_coordinate(i));
     if (cr != EQUAL) return cr;
   }
   return EQUAL;
 }

  CGAL_AOS3_EVENT_1K(I, insert);
  CGAL_AOS3_EVENT_1K(R, remove);

  CGAL_AOS3_EVENT_2K(I2, intersect);
  CGAL_AOS3_EVENT_2K(U2, unintersect);

  CGAL_AOS3_EVENT_2K(S2, swap);

  CGAL_AOS3_EVENT_3K(I3, intersect);
  CGAL_AOS3_EVENT_3K(U3, unintersect);
  
  CGAL_AOS3_EVENT_1H(AAR, process_aar, is_sphere_sphere, is_sphere_rule);
  CGAL_AOS3_EVENT_1H(RAR, process_rar, is_sphere_rule, is_sphere_rule);
  CGAL_AOS3_EVENT_1H(ARR, process_arr, is_sphere_rule, is_rule_rule);

  CGAL_AOS3_EVENT_2K1I(AAE, process_aae);

  Event_processor(CCS &cs, Traits &tr);

  template <class Oit>
  void gather_incident_faces(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k, 
			    Oit out) {
    std::vector<CGAL_AOS3_TYPENAME CCS::Face_handle> border_faces;
    CGAL_AOS3_TYPENAME CCS::Halfedge_handle h= cs_.a_halfedge(k);
    h=h->opposite();
    CGAL_assertion(!h->curve().is_inside());
    do {
      border_faces.push_back(h->face());
      border_faces.push_back(h->opposite()->face());
      h= cs_.next_edge_on_circle(h);
    } while (h != cs_.a_halfedge(k)->opposite());
    std::sort(border_faces.begin(), border_faces.end(), CGAL_AOS3_TYPENAME CCS::Handle_compare());
    border_faces.erase(std::unique(border_faces.begin(), border_faces.end(),
				   CGAL_AOS3_TYPENAME CCS::Handle_equal()), border_faces.end());
    std::copy(border_faces.begin(), border_faces.end(), out);
  }


  void delete_free_event() {
    cs_.visitor().delete_free_event();
  }

  /*void insert(Sphere_3_key k);
  void remove(Sphere_3_key k);
  void intersect(Sphere_3_key k, Sphere_3_key l);
  void unintersect(Sphere_3_key k, Sphere_3_key l);

  void intersect(Sphere_3_key k, Sphere_3_key l,  Sphere_3_key m);
  void unintersect(Sphere_3_key k, Sphere_3_key l, Sphere_3_key m);*/

private:

  void check_degeneracy();

  void handle_degeneracy();

  void clear_event(Halfedge_handle h);
  CCS &cs_;
  Traits tr_;
};


CGAL_AOS3_END_INTERNAL_NAMESPACE


#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Event_processor_impl.h>
#endif

#endif
