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

#ifndef CGAL_KINETIC_IO_INTERNAL_QT_SIMULATOR_CORE_H
#define CGAL_KINETIC_IO_INTERNAL_QT_SIMULATOR_CORE_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/IO/internal/GUI_base.h>
#include <CGAL/Kinetic/Multi_listener.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <vector>
#include <qobject.h>

namespace CGAL
{
  namespace Kinetic
  {
    namespace internal
    {
      class Qt_core: public QObject, public Non_ref_counted<Qt_core>
      {
	typedef Qt_core This;
	Q_OBJECT
	public:
	//typedef int Timer_id;

	Qt_core();

	CGAL_KINETIC_LISTENERNT1(LAST_BUTTON_PRESSED)
      public:
      
	/*class Listener
	  {
	  public:
	  typedef Qt_core* Notifier_handle;
	  Listener(Qt_core *c): h_(c) {
	  c->set_listener(this);
	  }
	  virtual ~Listener() {
	  h_->set_listener(NULL);
	  }
	  typedef enum {LAST_BUTTON_PRESSED}
	  Notification_type;
	  virtual void new_notification(Notification_type tp)=0;
	  Qt_core *notifier() {
	  return h_;
	  }
	  protected:
	  Qt_core *h_;
	  };*/

	enum Button {RUN, STOP, RUN_TO, RUN_THROUGH, REVERSE, PAUSE, FASTER, SLOWER};
	typedef enum Button Button;

	Button last_button_pressed() const
	{
	  return mode_;
	}
      protected:
	//Listener *playable_;
	Button mode_;
      private:
	//friend class Listener;
	/*void set_listener(Listener *p) {
	  playable_= p;
	  }
	  const Listener *listener() const
	  {
	  return playable_;
	  }*/
      public slots:                     //functions
	
	void play_button();
	void pause_button();
	void stop_button();
	void play_to_button();
	void play_through_button();
	void reverse_button();
	void faster_button();
	void slower_button();
      };

      template <class Base>
      class Qt_core_listener: public Qt_core::Listener
      {
	typedef typename Qt_core::Listener IF;
	typedef Qt_core BH;
      public:
	Qt_core_listener(typename Base::Handle &t): t_(t) {
	}
	virtual void new_notification(typename IF::Notification_type nt) {
	  if (nt == IF::LAST_BUTTON_PRESSED) {
	    if (notifier()->last_button_pressed() == BH::RUN) {
	      t_->set_mode(Base::RUNNING);
	    }
	    else if (notifier()->last_button_pressed() == BH::PAUSE) {
	      t_->set_mode(Base::PAUSED);
	    }
	    else if (notifier()->last_button_pressed() == BH::STOP) {
	      t_->set_mode(Base::STOPPED);
	    }
	    else if (notifier()->last_button_pressed() == BH::RUN_TO) {
	      t_->set_mode(Base::RUNNING_TO_EVENT);
	    }
	    else if (notifier()->last_button_pressed() == BH::RUN_THROUGH) {
	      t_->set_mode(Base::RUNNING_THROUGH_EVENT);
	    }
	    else if (notifier()->last_button_pressed() == BH::REVERSE) {
	      if (t_->simulator()->direction_of_time()==CGAL::POSITIVE) {
		t_->simulator()->set_direction_of_time(CGAL::NEGATIVE);
	      }
	      else {
		t_->simulator()->set_direction_of_time(CGAL::POSITIVE);
	      }
	    }
	    else if (notifier()->last_button_pressed() == BH::FASTER) {
	      t_->set_speed(t_->speed()+.25);
	    }
	    else if (notifier()->last_button_pressed() == BH::SLOWER) {
	      t_->set_speed(t_->speed()-.25);
	    }
	  }
	}
	virtual ~Qt_core_listener() {
	}
      protected:
	typename Base::Handle t_;
      };

    }
  }
}
#endif
