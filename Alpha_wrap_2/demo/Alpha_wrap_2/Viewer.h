// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//                 Cedric Portaneri,
//                 Mael Rouxel-Labbé
//                 Andreas Fabri
//                 Michael Hemmer
//
#ifndef CGAL_ALPHA_WRAP_2_DEMO_VIEWER_H
#define CGAL_ALPHA_WRAP_2_DEMO_VIEWER_H

#include "scene.h"

#include <CGAL/Qt/qglviewer.h>

#include <QOpenGLWidget>
#include <QPaintEvent>

class Viewer
  : public CGAL::QGLViewer
{
private:
  Scene* m_scene;

  // camera
  double m_scale;
  double m_center_x, m_center_y;

  // mouse
  QPoint m_mouse_click, m_mouse_move;

  bool is_drawing = false;

public:
  Viewer(QWidget *parent);

  void set_scene(Scene* pScene) { m_scene = pScene; }

  void set_camera(const double x, const double y, const double s)
  {
    m_center_x = x;
    m_center_y = y;
    m_scale = s;
  }

protected:
  // GL
  void paintGL();
  void initializeGL();
  void resizeGL(int width, int height);

  // mouse
  void wheelEvent(QWheelEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);

  void keyPressEvent(QKeyEvent *event);
  void keyReleaseEvent(QKeyEvent *event);

  void sample_mouse_path(const QPoint& point, bool new_cmp, bool is_closed);
  void move_camera(const QPoint& p0, const QPoint& p1);
  void convert_to_world_space(const QPoint& point, double &x, double &y);
};

#endif // CGAL_ALPHA_WRAP_2_DEMO_VIEWER_H
