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

#ifndef CGAL_KINETIC_IO_QT_SIMULATOR_3_H_
#define CGAL_KINETIC_IO_QT_SIMULATOR_3_H_
#include <CGAL/Kinetic/basic.h>
#include "SoQt_handle.h"
#include <CGAL/Kinetic/IO/internal/GUI_base.h>
#include "SoQt_examiner_viewer.h"
#include <CGAL/Kinetic/IO/internal/Qt_timer.h>
#include <CGAL/Kinetic/Listener.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <Inventor/Qt/SoQt.h>
#include <Inventor/nodes/SoSeparator.h>

class SoSeparator;

namespace CGAL { namespace Kinetic {;

//! This provides a GUI in 3D using Coin.
/*!  If you want to draw something, extend
  SoQt_widget_3::Listener. SoQt_widget_3::Listener, in addition to the fields of
  CGAL::Listener, has an extra field root() which provides a
  SoSeperator node to act as the root of any scene graph nodes you
  wish to use.

  The GUI uses the SoQt widget set so Qt must be installed (it is
  needed for 2D also).

  An example using this GUI and SoQt_moving_points_3 and
  SoQt_moving_weighted_points_3 can be found in \example 3d_gui.cc.
*/
template <class Simulator_t>
class SoQt_widget_3:
  public Ref_counted<SoQt_widget_3<Simulator_t> >
{
protected:
  typedef SoQt_widget_3<Simulator_t> This;
  typedef Gui_base<Simulator_t,  internal::Qt_timer> Graphical_base;
  typedef typename Simulator_t::Time Time;
  typedef typename internal::Qt_core_listener<Graphical_base> Window_listener;
public:

  typedef Simulator_t Simulator;

  //! construct things
  SoQt_widget_3(int argc, char *argv[], typename Simulator::Handle sh):
    base_(new Graphical_base(sh)),
    window_l_(base_){
    CGAL_KINETIC_INIT_LISTEN(Graphical_base, base_);
    main_window_= SoQt::init(argc, argv, argv[0]);
    viewer_= new SoQt_examiner_viewer(main_window_);
    SoQt::show(main_window_);
    window_l_.set_notifier(viewer_->button_handler());
  }

  virtual ~SoQt_widget_3(){}

  //! start the gui
  int begin_event_loop() {
    update_coordinates();
    SoQt::mainLoop();
    return 0;
  }

  //! Return a (reference counted) pointer to the simulator.
  typename Simulator::Handle& simulator() {
    return base_->simulator();
  }
  //! Return a const (reference counted) pointer to the simulator
  typename Simulator::Handle simulator() const
  {
    return base_->simulator();
  }

  //! Get the current time as a double.
  double current_time() const
  {
    return base_->current_time();
  }

  SoQt_examiner_viewer* soqt_examiner_viewer_pointer() const  {
    return  viewer_;
  }

  /*class Listener_core
  {
  public:
    typedef typename This::Handle Notifier_handle;
    typedef enum {CURRENT_TIME}
      Notification_type;

    SoSeparator* root() {
      return parent_.get();
    }

    Listener_core(){}
  private:
    friend class SoQt_widget_3<Simulator_t>;
    void set_root(SoSeparator* p) {
      parent_=SoQt_handle<SoSeparator>(p);
    }
    SoQt_handle<SoSeparator> parent_;
  };
  */
  //! Extend this object to listen for events.
  /*!  If you create an instance of this listener, you will
    automatically be subscribed.
  */
  //typedef Multi_listener<Listener_core> Listener;
  //friend class Multi_listener<Listener_core>;

  struct Listener_core{						
    typedef typename This::Handle Notifier_handle;		
    typedef enum {CURRENT_TIME} Notification_type;
    SoSeparator* root() const {
      return parent_.get();
    }
  private:
    friend class SoQt_widget_3<Simulator_t>;
    void set_root(SoSeparator* p) {
      parent_=SoQt_handle<SoSeparator>(p);
    }
    SoQt_handle<SoSeparator> parent_;
  };								
public:								
 
  typedef Multi_listener_base<Listener_core> Listener;
  friend class Multi_listener_base<Listener_core> ;
private:							
 void new_listener(Listener *sk) {				
   listeners_.push_back(sk);
   SoSeparator* sep= new SoSeparator;
   viewer_->new_subgraph(sep);
   sk->set_root(sep);
 }								
 void delete_listener(Listener *kds) {				
   for (unsigned int i=0; i< listeners_.size(); ++i){		
     if (listeners_[i] == kds) {				
       std::swap(listeners_[i], listeners_.back());		
       listeners_.pop_back();	
       viewer_->delete_subgraph(kds->root());
       return;							
     }								
   }								
 }								
 std::vector<Listener*> listeners_;

private:

  CGAL_KINETIC_LISTEN1(Graphical_base, CURRENT_TIME, update_coordinates());

  void update_coordinates() {
    CGAL_KINETIC_MULTINOTIFY(CURRENT_TIME);
  }
protected:
  typename Graphical_base::Handle base_;
  QWidget *main_window_;
  SoQt_examiner_viewer* viewer_;
  Window_listener window_l_;
};

} } //namespace CGAL::Kinetic;
#endif                                            // qt
