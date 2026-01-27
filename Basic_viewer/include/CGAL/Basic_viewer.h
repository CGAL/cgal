// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_BASIC_VIEWER_H
#define CGAL_BASIC_VIEWER_H

#include <CGAL/license/GraphicsView.h>
#include <CGAL/Graphics_scene.h>

// compatibility
#if defined(CGAL_USE_BASIC_VIEWER) && !defined(CGAL_USE_BASIC_VIEWER_QT)
#define CGAL_USE_BASIC_VIEWER_QT 1
#endif

#if defined(CGAL_USE_BASIC_VIEWER_QT) && !defined(CGAL_USE_BASIC_VIEWER)
#define CGAL_USE_BASIC_VIEWER 1
#endif

#if defined(CGAL_USE_BASIC_VIEWER_QT)
#include <CGAL/Qt/Basic_viewer.h>
// #elif defined(CGAL_USE_BASIC_VIEWER_GLFW)
// #include <CGAL/GLFW/Basic_viewer.h>
#else
namespace CGAL
{
  inline
  void draw_graphics_scene(const Graphics_scene&,
                           const char* ="CGAL Basic Viewer")
  {
    std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
  }
} // End namespace CGAL
#endif

#endif // CGAL_BASIC_VIEWER_H
