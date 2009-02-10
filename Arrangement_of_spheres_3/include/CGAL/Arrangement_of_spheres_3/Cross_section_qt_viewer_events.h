#ifndef CGAL_AOS3_CROSS_SECTION_QT_VIEWER_EVENTS_H
#define CGAL_AOS3_CROSS_SECTION_QT_VIEWER_EVENTS_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/IO/Qt_examiner_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Rational_cross_section.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arrangement_of_spheres_3/Event_processor.h>
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Cross_section_qt_viewer_events {
  CGAL_AOS3_TRAITS;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CS;
  typedef Rational_cross_section CGAL_AOS3_TARG RCS;
  typedef CGAL::Simple_cartesian<double> K;
public:
  typedef CGAL_AOS3_TYPENAME Traits::FT NT;
  Cross_section_qt_viewer_events(Traits tr,  const CS &cs,
				 CGAL::Layer l,
				 CGAL::Color col= CGAL::ORANGE): tr_(tr),
							       cs_(cs),
							       rcs_(cs, tr_),
							       layer_(l), 
							       color_(col)
  {}

  void operator()(NT z, Qt_examiner_viewer_2 *qtv) {
    rcs_.set_z(z);
    //t_.set_temp_sphere(T::Sphere_3(T::Point_3(0,0,z), 0));
    *qtv << layer_;
    QMutexLocker lock(qtv->lock());
    *qtv << Qt_examiner_viewer_2::Erase();
    *qtv << color_;
    qtv->set_updating_box(false);
  
  

    for (CGAL_AOS3_TYPENAME CS::Halfedge_const_iterator hit= cs_.halfedges_begin();
	 hit != cs_.halfedges_end(); ++hit){
      if (hit->curve().key().is_target()) continue;
      if (hit->event() != CGAL_AOS3_TYPENAME Traits::Event_key()) {
	if (hit->event() != cs_.visitor().simulator()->null_event()) {
	  *qtv << color_;
	} else {
	  *qtv << CGAL::Color(static_cast<unsigned char>(color_.red()*.5), 
			      static_cast<unsigned char>(color_.green()*.5),
			      static_cast<unsigned char>(color_.blue()*.5));
	}
	if (hit->curve().is_rule() && hit->curve().is_inside()){
	  
	  CGAL_AOS3_TYPENAME K::Point_2 t= display_point(hit->vertex()->point());
	  CGAL_AOS3_TYPENAME K::Point_2 s= display_point(hit->opposite()->vertex()->point());
	  
	  *qtv << CGAL_AOS3_TYPENAME K::Segment_2(t,s);
	} else if (hit->curve().is_arc() && hit->curve().is_inside()){
	  CGAL_AOS3_TYPENAME K::Point_2 t= display_point(hit->vertex()->point());
	  CGAL_AOS3_TYPENAME K::Point_2 s= display_point(hit->opposite()->vertex()->point());
	  CGAL_AOS3_TYPENAME Traits::Circle_2 c2= rcs_.circle(hit->curve().key());
	  qtv->new_circular_arc(c2, s, t);
	}
      }
    }
   
    for (CGAL_AOS3_TYPENAME CS::Visitor::Free_event_const_iterator it= cs_.visitor().free_events_begin();
	 it != cs_.visitor().free_events_end(); ++it) {
      const CGAL_AOS3_TYPENAME Event_processor CGAL_AOS3_TARG::Event_base &eb= 
	 cs_.visitor().simulator()->event<CGAL_AOS3_TYPENAME Event_processor CGAL_AOS3_TARG::Event_base>(*it);
      double d= CGAL::to_interval(cs_.visitor().simulator()->event_time(*it)).first;
      if (!eb.is_edge_event() && CGAL::abs(d-z)<1) {
	CGAL_AOS3_TYPENAME Traits::Event_point_3 ep= cs_.visitor().simulator()->event_time(*it);
	std::ostringstream oss;
	eb.write(oss);
	*qtv << color_;
	*qtv << CGAL_AOS3_TYPENAME K::Point_2(ep.approximate_coordinate(plane_coordinate(0)),
					      ep.approximate_coordinate(plane_coordinate(1)))
	     << oss.str().c_str();
							      
      }
    }
  
    qtv->set_is_dirty(true);
  }



private:

  CGAL_AOS3_TYPENAME K::Point_2 display_point(CGAL_AOS3_TYPENAME CS::Point pt) const {
    try {
      CGAL_AOS3_TYPENAME Traits::Sphere_point_3 sp= rcs_.sphere_point(pt);
      //std::cout << "Exact point is " << sp << std::endl;
      return CGAL_AOS3_TYPENAME K::Point_2(sp.approximate_coordinate(plane_coordinate(0)),
					   sp.approximate_coordinate(plane_coordinate(1)));
    } catch (CGAL_AOS3_TYPENAME Traits::Point_3 pt) {
      std::cout << "Point " << pt << " is no longer valid at " << rcs_.z() << std::endl;
      return CGAL_AOS3_TYPENAME K::Point_2(CGAL::to_double(pt[plane_coordinate(0).index()]),
					   CGAL::to_double(pt[plane_coordinate(1).index()]));
    }
  }

  Traits tr_;
  const CS &cs_;
  RCS rcs_;
  Layer layer_;
  CGAL::Color color_;
};

CGAL_AOS3_END_INTERNAL_NAMESPACE


#endif
