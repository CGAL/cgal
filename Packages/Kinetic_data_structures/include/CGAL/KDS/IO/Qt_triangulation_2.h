// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_IO_QT_TRIANGULATION_2_H
#define CGAL_KDS_IO_QT_TRIANGULATION_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/KDS/Ref_counted.h>
#include <CGAL/KDS/internal/tds_2_helpers.h>

CGAL_KDS_BEGIN_NAMESPACE

//! This class draws a Kinetic_Delaunay_2 triangulation to a Qt_gui_2.
/*!  The most recently created edges are colored green and the other
  edges are colored black. See kinetic_Delaunay_2.cc for a useage
  example. There are no public methods other than the constructor.
*/
template <class Kinetic_Delaunay, class Qt_gui, class Qt_mpt> 
class Qt_triangulation_2: public Ref_counted<Qt_triangulation_2<Kinetic_Delaunay, Qt_gui, Qt_mpt> > {
  typedef Qt_triangulation_2<Kinetic_Delaunay, Qt_gui, Qt_mpt> This;
  typedef internal::Triangulation_data_structure_helper_2<typename Kinetic_Delaunay::Triangulation::Triangulation_data_structure> TDS_helper;
public:
  //typedef Kinetic_Delaunay Kinetic_Delaunay;
  //typedef CGAL::Ref_counted_pointer<This> Pointer;

  Qt_triangulation_2(typename Kinetic_Delaunay::Pointer &kdel,
		     typename Qt_gui::Pointer &gui, 
		     typename Qt_mpt::Pointer &mps): mpt_(mps),
									 listener_(gui, this),
									 kdel_(kdel){
  }
 
protected:

  

  //! This class listens for redraw requests (PICTURE_IS_VALID becoming false)
  /*!
    It calls the draw method when it recieves a notification.
  */
  class Listener: public Qt_gui::Listener{
    typedef typename Qt_gui::Listener P;
  public:
    Listener(typename Qt_gui::Pointer &h, This *t): Qt_gui::Listener(h), t_(t){}
    virtual void new_notification(typename Qt_gui::Listener::Notification_type nt) {
      if (nt == Qt_gui::Listener::PICTURE_IS_VALID){
	t_->draw(*P::widget(), P::notifier()->current_time());
      }
    }
  protected:
    This *t_;
  };
  friend class Listener;
 
  //! Draw the triangulation.
  void draw( CGAL::Qt_widget &w, double t) const {
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
	 fit != tri.finite_edges_end(); ++fit){
      Static_point p0= tri.geom_traits().static_object(fit->first->vertex((fit->second+1)%3)->point());
      Static_point p1= tri.geom_traits().static_object(fit->first->vertex((fit->second+2)%3)->point());
      Static_segment ss(p0, p1);
      if (!TDS_helper::get_undirected_edge_label(*fit)){
	w << CGAL::Color(125,125,125);
      } else if (kdel_->visitor().contains(*fit) || kdel_->visitor().contains(TDS_helper::mirror_edge(*fit))){
	w<< CGAL::Color(0,255,0);
      } else {
	w << CGAL::Color(0,0,0);
      }
      w << ss;
    }
  }

  typename Qt_mpt::Pointer mpt_;
  Listener listener_;
  typename Kinetic_Delaunay::Pointer kdel_;
};

CGAL_KDS_END_NAMESPACE
#endif
