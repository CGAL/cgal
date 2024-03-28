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
//                 Kumar Snehasish <kumar.snehasish@gmail.com>
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>
//
#ifndef VIEWER_H
#define VIEWER_H

#include "typedefs.h"
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Qt/Basic_viewer.h>

class Viewer : public CGAL::Basic_viewer
{
  Q_OBJECT

  typedef CGAL::Basic_viewer Base;

public:
  Viewer(QWidget* parent);
  void setScene(Scene* scene_, bool doredraw=true);
  void keyPressEvent(QKeyEvent *e);
  virtual QString helpString() const;

public Q_SLOTS:
  void sceneChanged();

private:
  CGAL::Graphics_scene_options<LCC,
                               Dart_const_descriptor,
                               Dart_const_descriptor,
                               Dart_const_descriptor,
                               Dart_const_descriptor> m_gs_options;
  CGAL::Graphics_scene m_graphic_buffer;
  Scene* scene;
  bool m_previous_scene_empty;
};

#endif
