// Copyright (c) 2002-2004,2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Radu Ursu

#include <CGAL/basic.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Qt_widget_Delaunay_triangulation_2.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_get_point.h>

#include <qapplication.h>
#include <qmainwindow.h>

typedef CGAL::Cartesian<double>             Rep;
typedef CGAL::Point_2<Rep>                  Point_2;
typedef CGAL::Delaunay_triangulation_2<Rep> Delaunay;

Delaunay dt;

class My_Layer : public CGAL::Qt_widget_layer{
  void draw(){
    *widget << dt;
  }
};

class My_Window : public QMainWindow {
  Q_OBJECT
public:
  My_Window(int x, int y){
    widget = new CGAL::Qt_widget(this, "CGAL Qt_widget");
    setCentralWidget(widget);
    resize(x,y);
    widget->attach(&get_point);
    widget->attach(&v);
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),
            this, SLOT(get_new_object(CGAL::Object)));
    widget->set_window(0, 600, 0, 600);
  };
private:	//members
  CGAL::Qt_widget_get_point<Rep> get_point;
  My_Layer v;
  CGAL::Qt_widget *widget;
private slots:
  void get_new_object(CGAL::Object obj)
  {
    Point_2 p;
    if (CGAL::assign(p, obj)) {
      dt.insert(p);
    }
    widget->redraw();
  }
}; //endclass

//  moc_source_file : layer.cpp
#include "layer.moc"

int main( int argc, char **argv )
{
    QApplication app( argc, argv );
    My_Window *W = new My_Window(600,600);
    app.setMainWidget( W );
    W->show();
    return app.exec();
}
