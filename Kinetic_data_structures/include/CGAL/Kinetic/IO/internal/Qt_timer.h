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

#ifndef CGAL_KINETIC_IO_INTERNAL_QT_TIMER_H
#define CGAL_KINETIC_IO_INTERNAL_QT_TIMER_H
#include <map>
#include <qtimer.h>
#include <CGAL/Kinetic/Listener.h>
#include <CGAL/Kinetic/Ref_counted.h>
namespace CGAL
{
  namespace Kinetic
  {
    namespace internal
    {
      class Qt_timer: public QObject, public Non_ref_counted<Qt_timer>
      {
	Q_OBJECT
    	CGAL_KINETIC_LISTENERNT1(TICKS)
     
      public:
	Qt_timer();
	
	int ticks() const
	{
	  return tick_;
	}
	void clear() {
	  //CGAL_precondition(id_!=-1);
	  if (id_!= -1) timer_.killTimer(id_);
	  id_=-1;
	}
	void run(double time_in_seconds);
      protected:
	QTimer timer_;
	int tick_;
	int id_;
	
	
	 private slots:
	 void timerDone();
      };
    }
  }
}
#endif
