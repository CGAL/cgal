// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sophie Fei Che <fei@cis.udel.edu>
//
// File Description : Demo of CGAL 3D Triangulation package

#include "MainWindow.h"
#include <QApplication>
#include <CGAL/Qt/init_ogl_context.h>

int main(int argc, char** argv)
{
  CGAL::Qt::init_ogl_context(2, 1);

  QApplication app(argc, argv);

  app.setOrganizationDomain("inria.fr");
  app.setOrganizationName("INRIA");
  app.setApplicationName("3D Triangulation Demo");

  // Import resources from libCGALQt (Qt6).
  CGAL_QT_INIT_RESOURCES;

  MainWindow mw;
  mw.show();

  return app.exec();
}
