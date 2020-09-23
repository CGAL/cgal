// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Pierre Alliez, Camille Wormser
//
//******************************************************************************
// File Description : demo of AABB tree on polyhedral edge and facet primitives
//
//******************************************************************************

#include "MainWindow.h"
#include <QApplication>
#include <CGAL/Qt/resources.h>
#include <QMimeData>


int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  app.setOrganizationDomain("inria.fr");
  app.setOrganizationName("INRIA");
  app.setApplicationName("AABB tree demo");
  //for windows
#if (QT_VERSION >= QT_VERSION_CHECK(5, 3, 0))
  app.setAttribute(Qt::AA_UseDesktopOpenGL);
#endif

  // Import resources from libCGALQt (Qt5).
  CGAL_QT_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();
  QStringList args = app.arguments();
  args.removeAt(0);

  if(!args.empty() && args[0] == "--use-meta")
  {
    mainWindow.setAddKeyFrameKeyboardModifiers(::Qt::MetaModifier);
    args.removeAt(0);
  }

  Q_FOREACH(QString filename, args)
    mainWindow.open(filename);

  return app.exec();
}

#  include "Scene.cpp"
#  include "Scene_moc.cpp"
#  include "benchmarks.cpp"
#  include "Viewer.cpp"
#  include "Viewer_moc.cpp"
#  include "MainWindow.cpp"
#  include "MainWindow_moc.cpp"

