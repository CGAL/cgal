// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#include "MainWindow.h"
#include "typedefs.h"
#include <QApplication>
#include <CGAL/Qt/resources.h>
#include <CGAL/assertions_behaviour.h>
#include <CGAL/Qt/init_ogl_context.h>

// Global random
CGAL::Random myrandom;

int main(int argc, char** argv)
{
  // std::cout<<"Size of dart: "<<sizeof(LCC::Dart)<<std::endl;
  CGAL::set_error_behaviour(CGAL::ABORT);

  CGAL::Qt::init_ogl_context(4,3);
  QApplication application(argc,argv);

  application.setOrganizationDomain("cgal.org");
  application.setOrganizationName("CNRS and LIRIS' Establishments");
  application.setApplicationName("3D Linear Cell Complex");

  // Import resources from libCGALQt6
  // See https://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_Qt_init_resources();// that function is in a DLL
  Q_INIT_RESOURCE(Linear_cell_complex_3);

  MainWindow mw;
  mw.show();

  return application.exec();
}
