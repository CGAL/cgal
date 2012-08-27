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

#ifndef CGAL_KINETIC_IO_INTERNAL_QT_SIMULATOR_2_GUI_H
#define CGAL_KINETIC_IO_INTERNAL_QT_SIMULATOR_2_GUI_H
//#include <qtimer.h>
#include <CGAL/Kinetic/basic.h>

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/Kinetic/IO/internal/Qt_core.h>
#include <CGAL/Kinetic/IO/internal/pixmaps.h>
#include <CGAL/Kinetic/IO/internal/Qt_widget_2_core.h>
#include <map>
#include <qmainwindow.h>
#include <set>

// I think I need these here explicitly for MOC to work
namespace CGAL
{
  namespace Kinetic
  {
    namespace internal
    {
      /*
	Usage
	Qt_simulator_window win(-10,10, -10,10);
	QApplication app(argc, argv);
	app.setMainWidget( &_win );
	win.show();
	win.setCaption("KDS");
	app.exec();
      */
      class Qt_window_2 : public ::QMainWindow
      {
	Q_OBJECT
      public:

	~Qt_window_2(){}

	Qt_window_2(int xmin, int xmax, int ymin, int ymax);

	typedef Qt_core Button_handler;

	Button_handler* button_handler() {
	  return &core_;
	}

	Qt_widget_2_core *widget() {
	  return widget_;
	}

	Qt_widget_2_core *widget() const
	{
	  return widget_;
	}

	/*void redraw() // not redraw_win
	  {
	  //std::cout << "External redraw.\n";
	  _widget->redraw();
	  }*/

      private:                          //members
	CGAL::Qt_widget_standard_toolbar *_std_toolbar;
	Qt_widget_2_core *widget_;
	Qt_core core_;
      };
    }
  }
}
#endif
