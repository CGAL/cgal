// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

namespace CGAL { namespace Kinetic {

//! This class draws a Kinetic_Delaunay_2 triangulation to a Qt_gui_2.
/*!  The most recently created edges are colored green and the other
  edges are colored black. See kinetic_Delaunay_2.cc for a useage
  example. There are no public methods other than the constructor.
*/
template <class KDel, class IK, class Qt_gui>
class Qt_triangulation_2: public Ref_counted<Qt_triangulation_2<KDel, IK, Qt_gui> >
{
  typedef Qt_triangulation_2<KDel, IK, Qt_gui> This;
  // in icc Qt_gui::Listener::This captures This so I need to change the name
  //typedef This Qt_tri;;
  typedef typename KDel::Triangulation Triangulation;
  typedef internal::Triangulation_data_structure_helper_2<typename Triangulation::Triangulation_data_structure> TDS_helper;

  typedef typename Triangulation::Edge Edge;

  // maybe icl wants the class definition before the useage. 
  CGAL_KINETIC_LISTEN1(Qt_gui, PICTURE_IS_VALID, draw())


public:
  //typedef Kinetic_Delaunay Kinetic_Delaunay;
  //typedef CGAL::Ref_counted_pointer<This> Pointer;

  Qt_triangulation_2(typename KDel::Handle kdel,
		     IK ik,
		     typename Qt_gui::Handle gui): 
						   ik_(ik),
						   kdel_(kdel) {
    CGAL_KINETIC_INIT_LISTEN(Qt_gui, gui);
  }

protected:

  //! This class listens for redraw requests (PICTURE_IS_VALID becoming false)
  /*!
    It calls the draw method when it recieves a notification.
  */

  template <class V>
  void set_color(const Triangulation  &, 
                 const Edge &e, CGAL::Qt_widget &w, const V &) const {
    if (!kdel_->has_event(e)) {
      w << CGAL::Color(125,125,125);
    } else if (kdel_->has_finite_event(e)){
      w << CGAL::Color(255,0,0);
    } else {
      w << CGAL::Color(0,0,0);
    }
  }

  typedef Delaunay_triangulation_recent_edges_visitor_2<typename KDel::Triangulation> REV;

  void set_color(const Triangulation  &tri, 
                 const Edge &e, CGAL::Qt_widget &w,
                 const REV& ) const {
    w << CGAL::LineWidth(2);
    if (!kdel_->has_event(e)) {
      w << CGAL::Color(125,125,125);
    } else if (kdel_->visitor().contains(e) || kdel_->visitor().contains(tri.mirror_edge(e))) {
      w<< CGAL::Color(0,255,0);
    } else if (kdel_->has_finite_event(e)){
      w << CGAL::Color(125,125,125);
    } else {
      w << CGAL::Color(0,0,0);
    }
  }
  void draw() const {
    draw(CGAL_KINETIC_NOTIFIER(Qt_gui)->widget(), CGAL_KINETIC_NOTIFIER(Qt_gui)->current_time());
  }
  //! Draw the triangulation.
  void draw( CGAL::Qt_widget &w, double t) const
  {
   
    
    const Triangulation  &tri= kdel_->triangulation();
    ik_.set_time(typename IK::Time(t));
    typedef typename IK::Static_kernel::Point_2 Static_point;
    typedef typename IK::Static_kernel::Segment_2 Static_segment;
    w << CGAL::LineWidth(1);
    // << CGAL::FillColor(CGAL::Color(0,0,0));
    if (tri.dimension() != 2) return;
    ik_.set_time(typename IK::NT(t));
    typename IK::Current_coordinates cc= ik_.current_coordinates_object();
    for (typename Triangulation::Finite_edges_iterator fit = tri.finite_edges_begin();
	 fit != tri.finite_edges_end(); ++fit) {
      Static_point p0= cc(fit->first->vertex((fit->second+1)%3)->point());
      Static_point p1= cc(fit->first->vertex((fit->second+2)%3)->point());
      Static_segment ss(p0, p1);
      set_color(tri, *fit, w, kdel_->visitor());
      w << ss;
    }
  }

  IK ik_;
  typename KDel::Handle kdel_;
};

} } //namespace CGAL::Kinetic
#endif
