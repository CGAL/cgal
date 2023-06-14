// Copyright (C) 2016 The Qt Company Ltd.
// SPDX-License-Identifier: LicenseRef-Qt-Commercial OR BSD-3-Clause

#ifndef MAINWIDGET_H
#define MAINWIDGET_H

#include <QOpenGLWidget>
#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector2D>
#include <QBasicTimer>

#include <qopenglwidget.h>

#include <memory>

#include "Camera.h"
#include "Common_defs.h"
#include "Shader_program.h"
#include "Sphere.h"
#include "World_coordinate_axes.h"


class MainWidget : public QOpenGLWidget, protected OpenGLFunctionsBase
{
  Q_OBJECT

public:
  using QOpenGLWidget::QOpenGLWidget;
  ~MainWidget();

protected:
  void set_mouse_button_pressed_flag(QMouseEvent* e, bool flag);
  void mousePressEvent(QMouseEvent *e) override;
  void mouseMoveEvent(QMouseEvent* e) override;
  void mouseReleaseEvent(QMouseEvent *e) override;
  void timerEvent(QTimerEvent *e) override;

  void initializeGL() override;
  void resizeGL(int w, int h) override;
  void paintGL() override;

    
  void init_camera();
  void init_geometry();
  
  void init_shader_programs();
  void init_sp_smooth();
  void init_sp_color_only();

private:
  // Objects in the scene
  std::unique_ptr<Sphere>           m_sphere;
  std::unique_ptr<World_coord_axes> m_world_coord_axes;

  // Shaders
  Shader_program  m_sp_smooth;
  Shader_program  m_sp_color_only;
  
  // Camera & controls
  Camera m_camera;
  bool m_left_mouse_button_down = false;
  bool m_middle_mouse_button_down = false;
  QVector2D m_last_mouse_pos;

  // Timer for continuous screen-updates
  QBasicTimer m_timer;
};

#endif // MAINWIDGET_H
