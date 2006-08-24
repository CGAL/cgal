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

#ifndef CGAL_KINETIC_CGAL_GRAPHICAL_SIMULATOR_BASE_H
#define CGAL_KINETIC_CGAL_GRAPHICAL_SIMULATOR_BASE_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/Listener.h>

CGAL_KINETIC_BEGIN_NAMESPACE
//! This is a base class for implementing GUIs for kinetic simulations
/*!  Use this if you need to implement your on GUI. See Qt_gui_2 for a
  simple example of using it.
*/
template <class Simulator_t, class Timer>
class Gui_base: public CGAL::Kinetic::Ref_counted<Gui_base<Simulator_t, Timer> >
{
protected:
  typedef Simulator_t Simulator;
  typedef Gui_base<Simulator, Timer> This;

  typedef typename Simulator::Time TT;

public:
  //typedef CGAL::Ref_counted_pointer<This> Pointer;

  //! The current mode of the simulation.
  /*!
    The modes are:
    - RUNNING
    - RUNNING_TO_EVENT: run until just before the next event
    - RUNNING_THROUGH_EVENT: run until just after the next event
    - PAUSED: suspend whatever was being done before
    - STOPPED: there really is no function to this other than conventionality

    This creates notifications whenever the current time changes to
    allow the layer which interacts with the graphics to draw things.
  */
  typedef enum Mode {RUNNING, RUNNING_TO_EVENT, RUNNING_THROUGH_EVENT, PAUSED, STOPPED}
    Mode;

  //! Initialize with a pointer to the simulator.
  Gui_base(typename Simulator::Handle sh): mode_(STOPPED), paused_mode_(STOPPED),
					    fps_(30), speed_log_(0),
					    dir_of_time_(1), timer_(new Timer()),
					    timer_callback_(timer_,const_cast<This*>(this)),
					    drawable_(NULL), processing_(false) {
    sim_= sh;
  }

  virtual ~Gui_base() {
    delete timer_;
  }

  //! Return the current mode.
  Mode mode() {
    return mode_;
  }

  //! Set the current mode and update things accordingly.
  void set_mode(Mode mode) {
    if (mode== PAUSED) {
      if (mode_ == PAUSED) {
	unpause();
      }
      else {
	pause();
      }
    } else if ((mode==RUNNING_TO_EVENT || mode==RUNNING_THROUGH_EVENT)
	       && sim_->empty()) {
      std::cout << "No more events.\n";
      return;
    }
    else {
      mode_= mode;
    }

    CGAL_KINETIC_LOG(LOG_SOME, "Mode changed to " << mode_string(this->mode()) << std::endl);
    timer_->clear();
    switch(this->mode()) {
    case RUNNING_TO_EVENT:
    case RUNNING_THROUGH_EVENT:
    case RUNNING:
      //std::cout << "Starting timer with time " << Parent::step_length_real_time() << std::endl;
      timer_->run(step_length_real_time());
      break;
    case STOPPED:
    case PAUSED:
      break;
    default:
      std::cerr << "Invalid case: " << this->mode() << std::endl;
      CGAL_assertion(0);
    }
  }

  class Listener_core
  {
  public:
    typedef typename This::Handle Notifier_handle;
    typedef enum {CURRENT_TIME}
      Notification_type;
  };
  //! The class to extend if you want to receive events.
  /*!  See CGAL::Listener to a description of how runtime
    notifications are handled.
  */
  typedef CGAL::Kinetic::Listener<Listener_core> Listener;
  friend class CGAL::Kinetic::Listener<Listener_core>;

  //! get the simulator
  typename Simulator::Handle& simulator() {
    return sim_;
  }

  //! get the simulator
  typename Simulator::Const_handle simulator()const
  {
    return sim_;
  }

  //! The current time as a double.
  double current_time() const
  {
    return to_double(sim_->current_time());
  }

  //! The speed in logrithmic units.
  double speed() const
  {
    return speed_log_;
  }
  //! Set the speed of the simulation
  void set_speed(double s) {
    speed_log_=s;
  }

protected:

  const char * mode_string(Mode m) const
  {
    switch (m) {
    case RUNNING:
      return "<Run>";
      break;
    case RUNNING_TO_EVENT:
      return  "<Run to>";
      break;
    case RUNNING_THROUGH_EVENT:
      return "<Run through>";
      break;
    case PAUSED:
      return "<Pause>";
      break;
    case STOPPED:
      return "<Stop>";
      break;
    default:
      return "<Unknown>";
    }
  }

  const Listener* listener() const
  {
    return drawable_;
  }

  void set_listener(Listener* d) {
    drawable_=d;
  }

  class Timer_listener: public Timer::Listener
  {
  public:
    Timer_listener(Timer *tm, This *t):Timer::Listener(tm),  t_(t) {
    }
    void new_notification(typename Timer::Listener::Notification_type) {
      t_->timer_rang();
    }
  protected:
    This *t_;
  };

  friend class Timer_listener;

  void timer_rang() {
    // do something here
    if (increment_time()) {
      timer_->run(step_length_real_time());
    }
  }

  TT compute_next_timestep() const
  {
    return compute_next_timestep(sim_->end_time());
  }

  TT compute_next_timestep(const TT& max_time) const
  {
    double ct= to_double(sim_->current_time());
    double nt= ct + step_length_kds_time();
    //Parent::set_current_time(typename Parent::Time(nt));
    TT t(nt);
    if (CGAL::compare(t, sim_->current_time()) == CGAL::SMALLER) {
      std::cerr << "Time failed to advance due to roundoff. This is bad.\n";
      return sim_->current_time();
    }
    //std::cout << "Time is " << t << std::endl;
    if (CGAL::compare(t, max_time) == CGAL::LARGER) return max_time;
    else return t;
  }

  double step_length_real_time() const
  {
    return 1/fps_;
  }

  double step_length_kds_time() const
  {
    return  std::exp(speed_log_)/fps_;
  }

  bool increment_time() {
    CGAL_precondition(!is_processing());
    set_is_processing(true);
    //typename Parent::Time t= do_timestep(Parent::_end_time);
    bool ret(true);                       // to satisfy compiler
    switch (mode_) {
    case RUNNING:
      {
	TT nt= compute_next_timestep();
	sim_->set_current_time(nt);
	if (sim_->current_time() == sim_->end_time()) {
	  std::cout << "Simulation reached end of time. "
		    << "Press \"step over\" to process the next final event (if any)." << std::endl;
	  ret= false;
	}
	else ret= true;
	break;
      }
    case RUNNING_TO_EVENT:
      {
	TT stop_time= sim_->next_event_time();
	TT t= compute_next_timestep(stop_time);
	if (t == stop_time) {
	  set_mode(STOPPED);
	  //sim_->set_current_event_number(sim_->current_event_number()+1);
	  ret= false; break;
	}
	else {
	  sim_->set_current_time(t);
	  ret= true; break;
	}
      }
    case RUNNING_THROUGH_EVENT:
      {
	TT stop_time= sim_->next_event_time();
	TT t= compute_next_timestep(stop_time);
	if (t == stop_time) {
	  set_mode(STOPPED);
	  sim_->set_current_event_number(sim_->current_event_number()+1);
	  ret= false; break;
	}
	else {
	  sim_->set_current_time(t);
	  ret= true; break;
	}
      }
    default:
      std::cerr << "Run callback in invalid mode." << std::endl;
    }
    if (drawable_ != NULL) drawable_->new_notification(Listener::CURRENT_TIME);

    set_is_processing(false);

    return ret;
  }

  void pause() {
    paused_mode_= mode_;
    mode_=PAUSED;
  }
  void unpause() {
    mode_=paused_mode_;
  }

  bool is_processing() const
  {
    return processing_;
  }
  void set_is_processing(bool tf) {
    processing_=tf;
  }

  Mode mode_;
  Mode paused_mode_;
  double fps_;
  double speed_log_;
  int dir_of_time_;
  typename Simulator::Handle sim_;
  Timer *timer_;
  Timer_listener timer_callback_;
  Listener *drawable_;
  bool processing_;
};
CGAL_KINETIC_END_NAMESPACE;
#endif
