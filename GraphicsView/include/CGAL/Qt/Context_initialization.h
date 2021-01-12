// Copyright (c) 2021  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_QT_CONTEXT_INITIALIZATION_H
#define CGAL_QT_CONTEXT_INITIALIZATION_H
namespace CGAL
{
void init_ogl_context(int major, int minor) {
  QSurfaceFormat fmt;
#ifdef Q_OS_MAC
  fmt.setVersion(4, 1);
#else
  fmt.setVersion(4, 3);
#endif
  fmt.setRenderableType(QSurfaceFormat::OpenGL);
  fmt.setProfile(QSurfaceFormat::CoreProfile);
  fmt.setOption(QSurfaceFormat::DebugContext);
  QSurfaceFormat::setDefaultFormat(fmt);

  //for windows
#if (QT_VERSION >= QT_VERSION_CHECK(5, 3, 0))
  QCoreApplication::setAttribute(::Qt::AA_UseDesktopOpenGL);
#endif

  //We set the locale to avoid any trouble with VTK
  setlocale(LC_ALL, "C");
  return i;
}

}

#endif // CGAL_QT_CONTEXT_INITIALIZATION_H
