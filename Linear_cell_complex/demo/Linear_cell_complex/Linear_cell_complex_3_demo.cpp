// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#include "MainWindow.h"
#include "typedefs.h"
#include <QApplication>
#include <CGAL/Qt/resources.h>
#include <CGAL/assertions_behaviour.h>

// Global random
CGAL::Random myrandom;

int main(int argc, char** argv)
{
  // std::cout<<"Size of dart: "<<sizeof(LCC::Dart)<<std::endl;
  CGAL::set_error_behaviour(CGAL::ABORT);

  QApplication application(argc,argv);
  
  application.setOrganizationDomain("cgal.org");
  application.setOrganizationName("CNRS and LIRIS' Establishments");
  application.setApplicationName("3D Linear Cell Complex");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  CGAL_Qt4_init_resources(); // that function is in a DLL
  Q_INIT_RESOURCE(Linear_cell_complex_3);
  MainWindow mw;
  mw.show();

  return application.exec();
}
