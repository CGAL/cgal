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

#ifndef CGAL_KINETIC_IO_INTERNAL_QT_WIDGET_2_CORE_H
#define CGAL_KINETIC_IO_INTERNAL_QT_WIDGET_2_CORE_H

#include <CGAL/Kinetic/basic.h>
#include <CGAL/IO/Qt_widget.h>
#include <qmainwindow.h>
#include <CGAL/Kinetic/Listener.h>
#include <CGAL/Kinetic/Ref_counted.h>

namespace CGAL
{
  namespace Kinetic
  {
    namespace internal
    {
      class Qt_widget_2_core: public ::CGAL::Qt_widget, public Non_ref_counted<Qt_widget_2_core>
      {
	Q_OBJECT
      private:			
	struct Listener_core{						
	  typedef  This Notifier;		
	  typedef enum {PICTURE_IS_CURRENT} Notification_type;		
	};								
      public:							
	typedef CGAL::Kinetic::Listener_base<Listener_core> Listener;	
	friend class CGAL::Kinetic::Listener_base<Listener_core>;	
      private:							
	void set_listener(Listener *sk) {				
	  listener_=sk;						
	}								
	Listener* listener() {return listener_.get();}		
	Listener::Handle listener_;
      public:

	/*class Listener
	{
	public:
	  Listener(Qt_widget_2_core *widget): widget_(widget) {
	    CGAL_precondition(widget!= NULL);
	    widget_->set_listener(this);
	  }
	  virtual ~Listener() {
	    // could check first
	    widget_->set_listener(NULL);
	  }
	  typedef enum {PICTURE_IS_CURRENT}
	  Notification_type;
	  virtual void new_notification(Notification_type) {
	    //CGAL_ERROR( "draw not implemented.\n");
	    std::cerr << "Drawing but nothing is to be drawn.\n";
	  }
	  Qt_widget_2_core *widget(){return widget_;}
	protected:
	  Qt_widget_2_core *widget_;
	  };*/

	Qt_widget_2_core(QMainWindow *parent);

	//! do not call, this is for Qt use.
	void redraw() ;

	bool picture_is_current() const
	{
	  return is_drawn_;
	}
	void set_picture_is_current(bool tf) {
	  if (tf==false) redraw();
	}
      protected:
	bool is_drawn_;

      };
    }
  }
}
#endif
