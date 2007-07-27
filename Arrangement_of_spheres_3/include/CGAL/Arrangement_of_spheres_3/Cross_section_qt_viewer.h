#ifndef CGAL_AOS3_CROSS_SECTION_QT_VIEWER_H
#define CGAL_AOS3_CROSS_SECTION_QT_VIEWER_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/IO/Qt_examiner_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Rational_cross_section.h>
#include <CGAL/Simple_cartesian.h>
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Cross_section_qt_viewer {
  CGAL_AOS3_TRAITS;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CS;
  typedef Rational_cross_section CGAL_AOS3_TARG RCS;
  typedef CGAL::Simple_cartesian<double> K;
public:
  typedef CGAL_AOS3_TYPENAME Traits::FT NT;
  Cross_section_qt_viewer(Traits tr,  const CS &cs): tr_(tr),
										cs_(cs),
										rcs_(cs, tr_){}

  void operator()(NT z, Qt_examiner_viewer_2 *qtv) {
    rcs_.set_z(z);
    //t_.set_temp_sphere(T::Sphere_3(T::Point_3(0,0,z), 0));
    
    *qtv << CGAL::RED;
    qtv->set_updating_box(true);
    //T::Intersect_with_sweep is=t_.sphere_intersects_rule(z);
    
    
    /*for (CGAL_AOS3_TYPENAME Traits::Sphere_key_iterator sit= t_.sphere_keys_begin(); 
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

    for (CGAL_AOS3_TYPENAME CS::Halfedge_const_iterator hit= cs_.halfedges_begin();
	 hit != cs_.halfedges_end(); ++hit){
      if (hit->curve().key().is_target()) continue;
      if (hit->curve().is_rule() && hit->curve().is_inside()){
	qtv->set_updating_box(false);
	std::cout << "Displaying rule " << hit->curve() << std::endl;
	CGAL_AOS3_TYPENAME K::Point_2 t= display_point(hit->vertex()->point());
	CGAL_AOS3_TYPENAME K::Point_2 s= display_point(hit->opposite()->vertex()->point());
   
	*qtv << CGAL::GRAY;
	*qtv << CGAL_AOS3_TYPENAME K::Segment_2(t,s);
      } else if (hit->curve().is_arc() && hit->curve().is_inside()){
	qtv->set_updating_box(true);
	std::cout << "Displaying arc " << hit->curve() << std::endl;
	CGAL_AOS3_TYPENAME K::Point_2 t= display_point(hit->vertex()->point());
	CGAL_AOS3_TYPENAME K::Point_2 s= display_point(hit->opposite()->vertex()->point());
	//DT::Circle_2 c= ;
	if (tr_.compare_sphere_center_c(hit->curve().key(), z,
					sweep_coordinate())== CGAL::LARGER) {
	  *qtv << CGAL::Color(150,50,50);
	} else {
	  *qtv << CGAL::Color(50,150,50);
	}

	CGAL_AOS3_TYPENAME Traits::Circle_2 c2= rcs_.circle(hit->curve().key());
	qtv->new_circular_arc(c2, s, t);
      
      }
    }
   
    for (CGAL_AOS3_TYPENAME CS::Vertex_const_iterator hit= cs_.vertices_begin();
	 hit != cs_.vertices_end(); ++hit){
      if (!cs_.is_in_slice(hit)) continue;
      CGAL_AOS3_TYPENAME K::Point_2 p= display_point(hit->point());
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
    /*
     *qtv << CGAL::Color(255, 155, 155);
     for (CGAL_AOS3_TYPENAME Intersections_2::const_iterator it= intersections_2_.begin();
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
    */
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
      return CGAL_AOS3_TYPENAME K::Point_2(pt[plane_coordinate(0).index()],
					   pt[plane_coordinate(1).index()]);
    }
  }




  Traits tr_;
  const CS &cs_;
  RCS rcs_;
};

CGAL_AOS3_END_INTERNAL_NAMESPACE


#endif
