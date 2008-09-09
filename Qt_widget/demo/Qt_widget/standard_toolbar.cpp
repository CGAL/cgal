// Copyright (c) 2002-2004,2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent Rineau and Radu Ursu

#include <CGAL/basic.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Qt_widget_Delaunay_triangulation_2.h>

#include <qapplication.h>
#include <qmainwindow.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Cartesian<double>             Rep;
typedef CGAL::Point_2<Rep>                  Point_2;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;

Delaunay dt;

class My_window : public QMainWindow{
  Q_OBJECT
public:
  My_window(int x, int y)
  {
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);
    resize(x,y);
    widget->show();
    widget->set_window(0, x, 0, y);

    CGAL::Random_points_in_disc_2<Point_2> g(500);
    for(int count=0; count<100; count++) {
      dt.insert(*g++);
    }

    //How to attach the standard toolbar
    std_toolbar = new CGAL::Qt_widget_standard_toolbar(widget, this,
						       "Standard Toolbar");

    connect(widget, SIGNAL(redraw_on_back()),
	    this, SLOT(redraw_win()) );
  }

private slots:	//functions
  void redraw_win()
  {
    *widget << dt;
  }

private:	//members
  CGAL::Qt_widget *widget;
  CGAL::Qt_widget_standard_toolbar *std_toolbar;
};

// moc_source_file: standard_toolbar.cpp
#include "standard_toolbar.moc"

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_window W(600,600);
    app.setMainWidget( &W );
    W.show();
    W.setCaption("Using the Standard Toolbar");
    return app.exec();
}
