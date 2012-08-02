#ifndef QT_KINETIC_DELAUNAY_SS_2_H
#define QT_KINETIC_DELAUNAY_SS_2_H

#include <CGAL/IO/Qt_widget.h>
//#include <CGAL/Kinetic/IO/Qt_widget_2.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/internal/tds_2_helpers.h>
#include <CGAL/Kinetic/Listener.h>

namespace CGAL { namespace Kinetic {

//! This class draws a Kinetic_Delaunay_2 triangulation to a Qt_gui_2.
/*!  The most recently created edges are colored green and the other
  edges are colored black. See kinetic_Delaunay_2.cc for a useage
  example. There are no public methods other than the constructor.
*/
template <class Kinetic_Delaunay, class Qt_gui, class Qt_mpt>
class Qt_Delaunay_stable_subset_2: public Ref_counted<Qt_Delaunay_stable_subset_2<Kinetic_Delaunay, Qt_gui, Qt_mpt> >
{
  typedef Qt_Delaunay_stable_subset_2<Kinetic_Delaunay, Qt_gui, Qt_mpt> This;
  typedef Qt_Delaunay_stable_subset_2<Kinetic_Delaunay, Qt_gui, Qt_mpt> Qt_del;
  typedef internal::Triangulation_data_structure_helper_2<typename Kinetic_Delaunay::Triangulation::Triangulation_data_structure> TDS_helper;
public:
  //typedef Kinetic_Delaunay Kinetic_Delaunay;
  //typedef CGAL::Ref_counted_handle<This> Pointer;

  Qt_Delaunay_stable_subset_2(typename Qt_gui::Handle &gui,
			      typename Qt_mpt::Handle &mps,
			      typename Kinetic_Delaunay::Handle &kdel,
			      double threshold): mpt_(mps),
						 kdel_(kdel),
						 threshold_(threshold) {
    CGAL_KINETIC_INIT_LISTEN(Qt_gui, gui);
  }

protected:
  typedef typename Kinetic_Delaunay::Triangulation::Geom_traits::Static_kernel::Point_2 Static_point;
  typedef typename Kinetic_Delaunay::Triangulation::Geom_traits::Static_kernel::Vector_2 Static_vector;
  typedef typename Kinetic_Delaunay::Triangulation::Geom_traits::Static_kernel::Segment_2 Static_segment;
  typedef typename Kinetic_Delaunay::Triangulation::Geom_traits::Static_kernel::RT NT;

  double angle(const Static_point &a, const Static_point &b, const Static_point &c) const
  {
    Static_vector va(a-b);
    Static_vector vb(c-b);
    double dot= CGAL::to_double(va*vb);
    double ma= std::sqrt(CGAL::to_double(va*va));
    double mb= std::sqrt(CGAL::to_double(vb*vb));
    double ac= dot/(ma*mb);
    //std::cout << ac << " " << ma << " " << mb << std::endl;
    //std::cout << "Angle for " << a << ", " << b << ", " << c  << " is " << std::acos(ac) << std::endl;
    return std::acos(ac);
  }

  CGAL_KINETIC_LISTEN1(Qt_gui, PICTURE_IS_VALID, draw())
  
        typedef typename Qt_gui::Listener QTL;

  //! This class listens for redraw requests (PICTURE_IS_VALID becoming false)
  /*!
    It calls the draw method when it recieves a notification.
  
  class Listener: public QTL
  {
    typedef QTL P;
  public:
    Listener(typename Qt_gui::Handle &h, Qt_del *t): P(h), t_(t){}
    virtual void new_notification(typename QTL::Notification_type nt) {
      if (nt == QTL::PICTURE_IS_VALID) {
        t_->draw(*P::widget(), P::notifier()->current_time());
      }
    }
  protected:
    Qt_del *t_;
  };
  friend class Listener;*/

  void draw() const {
    draw(CGAL_KINETIC_NOTIFIER(Qt_gui)->widget(), CGAL_KINETIC_NOTIFIER(Qt_gui)->current_time());
  }

  //! Draw the triangulation.
  void draw( CGAL::Qt_widget &w, double t) const
  {
    //std::cout << "Drawing del\n";
    typedef typename Kinetic_Delaunay::Triangulation Del;
    const Del  &tri= kdel_->triangulation(typename Del::Geom_traits::Time(t));
    //tri.geom_traits().set_time(typename Del::Geom_traits::Time(t));
    w << CGAL::LineWidth(1);
    // << CGAL::FillColor(CGAL::Color(0,0,0));
    if (tri.dimension() != 2) return;
    for (typename Del::Finite_edges_iterator fit = tri.finite_edges_begin();
	 fit != tri.finite_edges_end(); ++fit) {
      if (fit->first->vertex((fit->second+1)%3)->point().is_valid()
	  && fit->first->vertex((fit->second+2)%3)->point().is_valid()
	  && fit->first->vertex(fit->second)->point().is_valid()
	  && fit->first->neighbor(fit->second)->vertex(tri.mirror_index(fit->first,fit->second))->point().is_valid()) {
	Static_point o= tri.geom_traits().current_coordinates_object()(fit->first->vertex((fit->second+1)%3)->point());
	Static_point d= tri.geom_traits().current_coordinates_object()(fit->first->vertex((fit->second+2)%3)->point());
	Static_point a= tri.geom_traits().current_coordinates_object()(fit->first->vertex(fit->second)->point());
	Static_point b= tri.geom_traits().current_coordinates_object()(fit->first->neighbor(fit->second)->vertex(tri.mirror_index(fit->first,fit->second))->point());

	double angle1= std::abs(angle(o,a, d));
	double angle2= std::abs(angle(o,b, d));

	Static_segment ss(o,d);

	if (angle1+angle2 < threshold_*3.1415) {
	  if (kdel_->visitor().contains(*fit) || kdel_->visitor().contains(tri.mirror_edge(*fit))) {
	    w<< CGAL::Color(255,0,0);
	  }
	  else {
	    w << CGAL::Color(0,0,0);
	  }
	  w << ss;
	} else {
	  Static_point o= tri.geom_traits().current_coordinates_object()(fit->first->vertex((fit->second+1)%3)->point());
	  Static_point d= tri.geom_traits().current_coordinates_object()(fit->first->vertex((fit->second+2)%3)->point());
	  Static_segment ss(o,d);
	  w << CGAL::Color(0,0,0);
	  w << ss;
	}
      }
      else {
	/*if (kdel_->visitor().contains(*fit) || kdel_->visitor().contains(TDS_helper::mirror_edge(*fit))){
	  w<< CGAL::Color(255,200,200);
	  } else {
	  w << CGAL::Color(200,200,200);
	  }*/
	//w << ss;
      }

      //Static_point p0= tri.geom_traits().current_coordinates_object()(fit->first->vertex((fit->second+1)%3)->point());
      //Static_point p1= tri.geom_traits().current_coordinates_object()(fit->first->vertex((fit->second+2)%3)->point());

    }
  }

  typename Qt_mpt::Handle mpt_;
  //Listener listener_;
  typename Kinetic_Delaunay::Handle kdel_;
  double threshold_;
};

} } //namespace CGAL::Kinetic
#endif
