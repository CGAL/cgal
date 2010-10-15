// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
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
#include <CGAL/IO/Qt_widget.h>

#include <qapplication.h>

typedef CGAL::Cartesian<int> Rep;
typedef CGAL::Point_2<Rep> Point_2;
typedef CGAL::Segment_2<Rep> Segment;

int main( int argc, char **argv )
{
  QApplication app( argc, argv );
  CGAL::Qt_widget *w = new CGAL::Qt_widget();
  app.setMainWidget( w );
  w->resize(600, 600);
  w->set_window(0, 600, 0, 600);
  w->show();
  w->lock();
  *w << CGAL::BackgroundColor(CGAL::ORANGE) << CGAL::RED;
  *w << Segment(Point_2(100,100), Point_2(400,400));
  w->unlock();
  return app.exec();
}
