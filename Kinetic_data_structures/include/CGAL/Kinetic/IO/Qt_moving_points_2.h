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

#ifndef CGAL_KINETIC_IO_MOVING_POINT_SET_2_H
#define CGAL_KINETIC_IO_MOVING_POINT_SET_2_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/Simulator_objects_listener.h>
#include <CGAL/Kinetic/Listener.h>
#include <sstream>

namespace CGAL { namespace Kinetic {

//! A graphical moving point set
/*!  This allows a set of moving points to be used in a
  CGAL::Kinetic::Qt_gui_2 GUI. It listens for PICTURE_IS_VALID
  notifications from the CGAL::Kinetic::Qt_gui_2 and DIRECTION_OF_TIME
  notifcations from the CGAL::Kinetic::Simulator and handles them
  accordingly.

*/
template <class TraitsC, class GUI>
class Qt_moving_points_2: public Ref_counted<Qt_moving_points_2<TraitsC, GUI> >
{
protected:
  typedef TraitsC Traits;
  //typedef Moving_object_set<Moving_point_t> Parent;
  typedef Qt_moving_points_2<Traits, GUI> This;

  typedef typename Traits::Simulator Simulator;
  // typedef typename GUI::Listener Gui_listener;
  //typedef typename Simulator::Listener Simulator_listener;
  //typedef CGAL::Kinetic::Simulator_objects_listener<Simulator_listener, This> Siml;
  //  friend class CGAL::Kinetic::Simulator_objects_listener<Simulator_listener, This>;
  typedef typename Traits::Static_kernel::FT NT;
  //class Guil;
  CGAL_KINETIC_LISTEN1(Simulator, DIRECTION_OF_TIME, reverse_time())
public:
  //typedef typename CGAL::Ref_counted_pointer<This> Pointer;
  //typedef typename CGAL::Const_ref_counted_pointer<This> Const_handle;

  //! Only outline drawing is currently supported
  typedef enum Draw_mode {FILLED, POINT, OUTLINE}
    Draw_mode;

  //! Defaults to outline drawing
  Qt_moving_points_2(typename GUI::Handle sim,
		     Traits tr): traits_(tr),
				 ik_(tr.instantaneous_kernel_object()),
				 cc_(ik_.current_coordinates_object()),
				 _mode(POINT), _radius(.2),
				 direction_of_time_(CGAL::POSITIVE),
				 //guil_(sim, const_cast<This*>(this)),
				 //siml_(tr.simulator_handle(),const_cast<This*>( this)),
				 rt_(tr.kinetic_kernel_object().reverse_time_object()) {
    CGAL_KINETIC_INIT_LISTEN(GUI, sim);
    CGAL_KINETIC_INIT_LISTEN(Simulator, tr.simulator_handle());
    ptsz_=10;
  };

  virtual ~Qt_moving_points_2(){}

  int point_size() const {
    return ptsz_;
  }
  
  void set_point_size(int p) {
    ptsz_ = p;
  }

  //! Change how things are drawn
  void set_draw_mode(Draw_mode mode) {
    _mode=mode;
  }

  //! change how big things are drawn
  void set_radius(double radius) {
    _radius=radius;
  }

  //! Set the CGAL::Sign field corresponding to the direction of time.
  /*!
    This can have the effect of reversing the motions of the points.
  */
  CGAL::Sign direction_of_time() const
  {
    return direction_of_time_;
  }


  unsigned int number_visible(NT xmin, NT xmax, NT ymin, NT ymax) const {
    typename Traits::Instantaneous_kernel::NT ntt(CGAL_KINETIC_NOTIFIER(GUI)->current_time());
    ik_.set_time(ntt);
    unsigned int ret=0;
    for (typename Traits::Active_points_2_table::Key_iterator
	   it= traits_.active_points_2_table_handle()->keys_begin();
	 it != traits_.active_points_2_table_handle()->keys_end(); ++it) {
      //std::cout << "drawing point " << *it  << "= " << ik_.to_static(*it) << std::endl;
      
      typename Traits::Static_kernel::Point_2 pt= cc_(*it);
      if (pt.x() >= xmin && pt.x() <= xmax
	  && pt.y() >= ymin && pt.y() <= ymax) {
	++ ret;
      }
    }
    return ret;
  }

  /*virtual void write(std::ostream &out) const {
    ik_.set_time(sim_->current_time_nt());
    for (typename Moving_point_table::key_iterator it= Moving_point_table::begin_keys(); it != Moving_point_table::end_keys(); ++it){
    out << *it;
    out << ": " << ik_.to_static(*it) << std::endl;
    }
    }*/

protected:
  CGAL_KINETIC_LISTEN1(GUI,PICTURE_IS_VALID, draw())
  /*class Guil: public Gui_listener
  {
    typedef Gui_listener P;
    typedef Qt_moving_points_2<Traits, GUI> QtMP;
  public:
    Guil(typename GUI::Handle &h, QtMP *t): Gui_listener(h), t_(t){}
    void new_notification(typename Gui_listener::Notification_type nt) {
      if (nt== Gui_listener::PICTURE_IS_VALID) {
	t_->draw();
      }
    }
  protected:
    QtMP *t_;
  };
  friend class Guil;*/

  void draw() const;

  void reverse_time();

  Traits traits_;
  typename Traits::Instantaneous_kernel ik_;
  typename Traits::Instantaneous_kernel::Current_coordinates cc_;
  Draw_mode _mode;
  double _radius;
  CGAL::Sign direction_of_time_;
  //Guil guil_;
  //Siml siml_;
  int ptsz_;
  typename Traits::Kinetic_kernel::Reverse_time rt_;
};

template <class T, class G>
void Qt_moving_points_2<T,G>::draw() const
{
  //std::cout << "Drawing mpt MPT\n";
  typedef typename Traits::Static_kernel::Circle_2 C;
  typedef typename Traits::Kinetic_kernel::Is_constant IC;
  IC ic = traits_.kinetic_kernel_object().is_constant_object();
  typename Traits::Simulator::NT ntt(CGAL_KINETIC_NOTIFIER(GUI)->current_time());
  ik_.set_time(ntt);

  CGAL::Qt_widget &w= CGAL_KINETIC_NOTIFIER(GUI)->widget();


  unsigned int numv= number_visible(w.x_min(), w.x_max(),
				    w.y_min(), w.y_max());
  if (numv < 200) {
    w << CGAL::PointSize(ptsz_) << CGAL::LineWidth(1);
  } else if (numv < 800) {
    w << CGAL::PointSize(std::max BOOST_PREVENT_MACRO_SUBSTITUTION(1,ptsz_/2)) << CGAL::LineWidth(1);
  } else if (numv < 2000) {
    w << CGAL::PointSize(std::max BOOST_PREVENT_MACRO_SUBSTITUTION(1,ptsz_/4)) << CGAL::LineWidth(1);
  } else {
    w << CGAL::PointSize(std::max BOOST_PREVENT_MACRO_SUBSTITUTION(ptsz_/8, 1)) << CGAL::LineWidth(1);
  }
  // typedef typename Traits::Static_kernel::Point_2 P2;
  //out << C(P2(0,0), 2) << C(P2(0,0), 1);
  //out << CGAL::BackgroundColor(CGAL::Color(125,125,125));
  
  for (typename Traits::Active_points_2_table::Key_iterator
	 it= traits_.active_points_2_table_handle()->keys_begin();
       it != traits_.active_points_2_table_handle()->keys_end(); ++it) {
    //std::cout << "drawing point " << *it  << "= " << ik_.to_static(*it) << std::endl;
    if (ic(traits_.active_points_2_table_handle()->at(*it))) {
      w << CGAL::FillColor(CGAL::Color(150,150,150));
      //out.setFillColor(CGAL::Color(0,255,0));
      w << CGAL::Color(150,150,150);
    } else {
      w << CGAL::FillColor(CGAL::Color(0,0,0));
      //out.setFillColor(CGAL::Color(0,255,0));
      w << CGAL::Color(0,0,0);
    }
    typename Traits::Static_kernel::Point_2 pt= cc_(*it);
    if (_mode== OUTLINE) {
      w << C(pt, _radius);
    }
    else if (_mode == POINT) {
      w << pt;
    }
    if (numv < 150) {
      std::ostringstream oss;
      oss << *it;
      w.get_painter().drawText(w.x_pixel(CGAL::to_double(pt.x()))+3,
				w.y_pixel(CGAL::to_double(pt.y()))-3,
				QString(oss.str().c_str()));
    }
  }
}


template <class T, class G>
void Qt_moving_points_2<T,G>::reverse_time()
{
  if (direction_of_time_ == CGAL::POSITIVE) direction_of_time_=CGAL::NEGATIVE;
  else direction_of_time_= CGAL::POSITIVE;

  traits_.active_points_2_table_handle()->set_is_editing(true);
  for (typename Traits::Active_points_2_table::Key_iterator kit
	 = traits_.active_points_2_table_handle()->keys_begin();
       kit != traits_.active_points_2_table_handle()->keys_end(); ++kit) {
    traits_.active_points_2_table_handle()->
      set(*kit, rt_(traits_.active_points_2_table_handle()->at(*kit)));
  }
  traits_.active_points_2_table_handle()->set_is_editing(false);
}


} } //namespace CGAL::Kinetic
#endif                                            // guard
