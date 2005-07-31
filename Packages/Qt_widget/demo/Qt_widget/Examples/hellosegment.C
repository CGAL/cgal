// Copyright (c) 1997-2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent Rineau
//                 Radu Ursu <rursu@sophia.inria.fr

#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*){
  std::cout << "Sorry, this demo needs QT..." << std::endl; return 0;}
#else
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
#endif
