// Copyright (C) 2016 The Qt Company Ltd.
// SPDX-License-Identifier: LicenseRef-Qt-Commercial OR BSD-3-Clause

#ifndef MAINWIDGET_H
#define MAINWIDGET_H

#include <memory>

#include <QOpenGLWidget>
#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector2D>
#include <QBasicTimer>

#include <qopenglwidget.h>
#include <qopenglfunctions_3_3_core.h>
#include <qopenglfunctions_4_5_core.h>
#include "Camera.h"

class Sphere;
using OpenGLFunctionsBase = QOpenGLFunctions_3_3_Core;


class MainWidget : public QOpenGLWidget, protected OpenGLFunctionsBase
{
    Q_OBJECT

public:
    using QOpenGLWidget::QOpenGLWidget;
    ~MainWidget();

protected:
    void mousePressEvent(QMouseEvent *e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent *e) override;
    void timerEvent(QTimerEvent *e) override;

    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;


    void add_shader(GLuint the_program, 
                    const char* shader_code, 
                    GLenum shader_type);
    void init_shader_program();
    
    void init_geometry();


private:

  std::unique_ptr<Sphere>  m_sphere;

  GLuint m_shader;
  GLuint m_uniform_mvp; // uniform location for MVP-matrix in the shader
  
  // camera & controls
  Camera m_camera;
  bool m_mouse_pressed = false;
  QVector2D m_last_mouse_pos;
  //QMatrix4x4 m_projection;


  QBasicTimer m_timer;
};

class Sphere : protected OpenGLFunctionsBase
{
public:
  Sphere(int num_slices, int num_stacks, float r);
  
  void draw();

private:
  GLuint m_vao, m_vbo, m_ibo, m_num_indices;
};

#endif // MAINWIDGET_H
