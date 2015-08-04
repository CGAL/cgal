// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://sccode@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/GSoC10-demoT3/demos/Triangulation_demo_3/ $
// $Id: r58170 2010-08-19 20:46:16 -0400 $
//
//
// Author(s)     : Sophie Fei Che <fei@cis.udel.edu>
//
// File Description : Demo of CGAL 3D Triangulation package

#include "MainWindow.h"
#include <QApplication>

int main(int argc, char** argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("inria.fr");
  app.setOrganizationName("INRIA");
  app.setApplicationName("3D Triangulation Demo");
  //for windows
#if (QT_VERSION >= QT_VERSION_CHECK(5, 3, 0))
  app.setAttribute(Qt::AA_UseDesktopOpenGL);
#endif

  MainWindow mw;
  mw.show();

  return app.exec();
}
