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

#ifndef CGAL_KDS_IO_INTERNAL_QT_SIMULATOR_CORE_H
#define CGAL_KDS_IO_INTERNAL_QT_SIMULATOR_CORE_H
#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/IO/internal/KDS_GUI_base.h>
#include <qobject.h>

namespace CGAL {
  namespace KDS {
    namespace internal {
      class Qt_core: public QObject{
	Q_OBJECT
      public:
	//typedef int Timer_id;

	Qt_core();

	class Listener {
	public:
	  typedef Qt_core* Notifier_pointer;
	  Listener(Qt_core *c): h_(c){
	    c->set_listener(this);
	  }
	  virtual ~Listener(){
	    h_->set_listener(NULL);
	  }
	  typedef enum {LAST_BUTTON_PRESSED} Notification_type;
	  virtual void new_notification(Notification_type tp)=0;
	  Qt_core *notifier() {
	    return h_;
	  }
	protected:
	  Qt_core *h_;
	};

	enum Button {RUN, STOP, RUN_TO, RUN_THROUGH, REVERSE, PAUSE, FASTER, SLOWER};
	typedef enum Button Button;

	Button last_button_pressed() const {
	  return mode_;
	}
      protected:
	Listener *playable_;
	Button mode_;
      private:
	friend class Listener;
	void set_listener(Listener *p){
	  playable_= p;
	}
	const Listener *listener() const {
	  return playable_;
	}
public slots:	//functions
   
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
      class Qt_core_listener: public Qt_core::Listener {
	typedef typename Qt_core::Listener IF;
	typedef Qt_core BH;
      public:
	Qt_core_listener(typename IF::Notifier_pointer w, typename Base::Pointer &t): IF(w),  t_(t){
	}
	virtual void new_notification(typename IF::Notification_type nt) {
	  if (nt == IF::LAST_BUTTON_PRESSED) {
	    if (notifier()->last_button_pressed() == BH::RUN){
	      t_->set_mode(Base::RUNNING);
	    } else if (notifier()->last_button_pressed() == BH::PAUSE){
	      t_->set_mode(Base::PAUSED);
	    } else if (notifier()->last_button_pressed() == BH::STOP){
	      t_->set_mode(Base::STOPPED);
	    } else if (notifier()->last_button_pressed() == BH::RUN_TO){
	      t_->set_mode(Base::RUNNING_TO_EVENT);
	    } else if (notifier()->last_button_pressed() == BH::RUN_THROUGH){
	      t_->set_mode(Base::RUNNING_THROUGH_EVENT);
	    } else if (notifier()->last_button_pressed() == BH::REVERSE){
	      if (t_->simulator()->direction_of_time()==CGAL::POSITIVE){
		t_->simulator()->set_direction_of_time(CGAL::NEGATIVE);
	      } else {
		t_->simulator()->set_direction_of_time(CGAL::POSITIVE);
	      }
	    } else if (notifier()->last_button_pressed() == BH::FASTER){
	      t_->set_speed(t_->speed()+.25);
	    } else if (notifier()->last_button_pressed() == BH::SLOWER){
	      t_->set_speed(t_->speed()-.25);
	    }
	  }
	}
	virtual ~Qt_core_listener(){
	}
      protected:
	typename Base::Pointer t_;
      };

    };
  };
};
#endif
