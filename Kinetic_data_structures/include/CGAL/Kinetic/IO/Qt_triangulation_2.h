// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_IO_QT_TRIANGULATION_2_H
#define CGAL_KINETIC_IO_QT_TRIANGULATION_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/internal/tds_2_helpers.h>
#include <CGAL/Kinetic/Delaunay_triangulation_recent_edges_visitor_2.h>

CGAL_KINETIC_BEGIN_NAMESPACE

//! This class draws a Kinetic_Delaunay_2 triangulation to a Qt_gui_2.
/*!  The most recently created edges are colored green and the other
  edges are colored black. See kinetic_Delaunay_2.cc for a useage
  example. There are no public methods other than the constructor.
*/
template <class Kinetic_Delaunay, class Qt_gui, class Qt_mpt>
class Qt_triangulation_2: public Ref_counted<Qt_triangulation_2<Kinetic_Delaunay, Qt_gui, Qt_mpt> >
{
  typedef Qt_triangulation_2<Kinetic_Delaunay, Qt_gui, Qt_mpt> This;
  // in icc Qt_gui::Listener::This captures This so I need to change the name
  //typedef This Qt_tri;;
  typedef internal::Triangulation_data_structure_helper_2<typename Kinetic_Delaunay::Triangulation::Triangulation_data_structure> TDS_helper;

  typedef typename TDS_helper::Edge Edge;

  // maybe icl wants the class definition before the useage. 
typedef typename Qt_gui::Listener QTL;
  class Listener: public QTL
  {
    typedef Qt_triangulation_2<Kinetic_Delaunay, Qt_gui, Qt_mpt> Container;
    typedef QTL P;
  public:
    Listener(typename Qt_gui::Handle &h, Container *t): P(h), t_(t){}
    virtual void new_notification(typename P::Notification_type nt) {
      if (nt == P::PICTURE_IS_VALID) {
	t_->draw(*P::widget(), P::notifier()->current_time());
      }
    }
  protected:
    Container *t_;
  };
  friend class Listener;

public:
  //typedef Kinetic_Delaunay Kinetic_Delaunay;
  //typedef CGAL::Ref_counted_pointer<This> Pointer;

  Qt_triangulation_2(typename Kinetic_Delaunay::Handle &kdel,
		     typename Qt_gui::Handle &gui,
		     typename Qt_mpt::Handle &mps): mpt_(mps),
						    listener_(gui, this),
						    kdel_(kdel) {
  }

protected:

  //! This class listens for redraw requests (PICTURE_IS_VALID becoming false)
  /*!
    It calls the draw method when it recieves a notification.
  */

  template <class V>
  void set_color(const Edge &e, CGAL::Qt_widget &w, const V &) const {
    if (!TDS_helper::get_undirected_edge_label(e)) {
      w << CGAL::Color(125,125,125);
    } else if (TDS_helper::get_undirected_edge_label(e) == kdel_->simulation_traits_object().simulator_handle()->null_event()){
      w << CGAL::Color(0,0,0);
    } else {
      w << CGAL::Color(255,0,0);
    }
  }

  typedef Delaunay_triangulation_recent_edges_visitor_2<typename Kinetic_Delaunay::Triangulation> REV;

  void set_color(const Edge &e, CGAL::Qt_widget &w,
		 const REV& ) const {
    w << CGAL::LineWidth(2);
    if (!TDS_helper::get_undirected_edge_label(e).is_valid()) {
      w << CGAL::Color(125,125,125);
    } else if (kdel_->visitor().contains(e) || kdel_->visitor().contains(TDS_helper::mirror_edge(e))) {
      w<< CGAL::Color(0,255,0);
    } else if (TDS_helper::get_undirected_edge_label(e) == kdel_->simulation_traits_object().simulator_handle()->null_event()){
      w << CGAL::Color(125,125,125);
    } else {
      w << CGAL::Color(0,0,0);
    }
  }

  //! Draw the triangulation.
  void draw( CGAL::Qt_widget &w, double t) const
  {
    //std::cout << "Drawing del\n";
    typedef typename Kinetic_Delaunay::Triangulation Del;

    const Del  &tri= kdel_->triangulation(typename Del::Geom_traits::Time(t));
    //tri.geom_traits().set_time(typename Del::Geom_traits::Time(t));
    typedef typename Del::Geom_traits::Static_kernel::Point_2 Static_point;
    typedef typename Del::Geom_traits::Static_kernel::Segment_2 Static_segment;
    w << CGAL::LineWidth(1);
    // << CGAL::FillColor(CGAL::Color(0,0,0));
    if (tri.dimension() != 2) return;
    for (typename Del::Finite_edges_iterator fit = tri.finite_edges_begin();
	 fit != tri.finite_edges_end(); ++fit) {
      Static_point p0= tri.geom_traits().static_object(fit->first->vertex((fit->second+1)%3)->point());
      Static_point p1= tri.geom_traits().static_object(fit->first->vertex((fit->second+2)%3)->point());
      Static_segment ss(p0, p1);
      set_color(*fit, w, kdel_->visitor());
      w << ss;
    }
  }

  typename Qt_mpt::Handle mpt_;
  Listener listener_;
  typename Kinetic_Delaunay::Handle kdel_;
};

CGAL_KINETIC_END_NAMESPACE
#endif
