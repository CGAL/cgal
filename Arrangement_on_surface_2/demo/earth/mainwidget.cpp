
#include "mainwidget.h"

#include <QMouseEvent>

#include <cmath>
#include <iostream>
#include <string>


MainWidget::~MainWidget()
{
  // Make sure the context is current when deleting the texture and the buffers.
  makeCurrent();
  doneCurrent();
}


void MainWidget::set_mouse_button_pressed_flag(QMouseEvent* e, bool flag)
{
  switch (e->button())
  {
  case Qt::LeftButton:
    m_left_mouse_button_down = flag;
    break;

  case Qt::MiddleButton:
    m_middle_mouse_button_down = flag;
    break;
  }
}
void MainWidget::mousePressEvent(QMouseEvent *e)
{
  set_mouse_button_pressed_flag(e, true);
  m_last_mouse_pos = QVector2D(e->position());
}
void MainWidget::mouseMoveEvent(QMouseEvent* e)
{
  auto current_mouse_pos = QVector2D(e->position());
  const auto diff = current_mouse_pos - m_last_mouse_pos;

  if (m_left_mouse_button_down)
  {  
    const float rotation_scale_factor = 0.1f;
    const float theta_around_x = rotation_scale_factor * diff.y();
    const float theta_around_y = rotation_scale_factor * diff.x();
    m_camera.rotate(theta_around_x, theta_around_y);
  }
  else if(m_middle_mouse_button_down)
  {
    const float zoom_scale_factor = 0.01f;
    const auto distance = zoom_scale_factor * diff.y();
    m_camera.move_forward(distance);
  }

  m_last_mouse_pos = current_mouse_pos;
}
void MainWidget::mouseReleaseEvent(QMouseEvent *e)
{
  set_mouse_button_pressed_flag(e, false);
}
void MainWidget::timerEvent(QTimerEvent *)
{
  update();
}



void MainWidget::initializeGL()
{
  initializeOpenGLFunctions();

  init_camera();
  init_geometry();
  init_shader_programs();

  glClearColor(0, 0, 0, 1);
  glEnable(GL_DEPTH_TEST);  // Enable depth buffer
  //glEnable(GL_CULL_FACE); // Enable back face culling

  // Use QBasicTimer because its faster than QTimer
  m_timer.start(12, this);
}



void MainWidget::init_camera()
{
  m_camera.set_pos(0, 0, 3);
  m_camera.rotate_around_x(-90);
}
void MainWidget::init_geometry()
{
  int num_slices, num_stacks;
  num_slices = num_stacks = 64;
  float r = 1;
  m_sphere = std::make_unique<Sphere>(num_slices, num_stacks, r);
  const float c = 0.8;
  m_sphere->set_color(c, c, c, 1);


  const float axes_length = 2;
  m_world_coord_axes = std::make_unique<World_coord_axes>(axes_length);
}
void MainWidget::init_shader_programs()
{
  init_sp_smooth();
  init_sp_color_only();
}
void MainWidget::init_sp_smooth()
{
  const char* vs = "shaders/smooth_vs.glsl";
  const char* fs = "shaders/smooth_fs.glsl";
  m_sp_smooth.init(vs, "", fs);

}
void MainWidget::init_sp_color_only()
{
  const char* vs = "shaders/color_only_vs.glsl";
  const char* fs = "shaders/color_only_fs.glsl";
  m_sp_color_only.init(vs, "", fs);
}



void MainWidget::resizeGL(int w, int h)
{
  // Reset projection
  qreal aspect = qreal(w) / qreal(h ? h : 1);
  const qreal z_near = 1.0, z_far = 100.0, fov = 45.0;
  m_camera.perspective(fov, aspect, z_near, z_far);
}
void MainWidget::paintGL()
{
  QMatrix4x4 model;
  const auto view = m_camera.get_view_matrix();
  const auto projection = m_camera.get_projection_matrix();
  const auto mvp = projection * view * model;

  // Clear color and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // SPHERE
  {
    auto& sp = m_sp_smooth;
    sp.use();
    sp.set_uniform("u_mvp", mvp);
    sp.set_uniform("u_color", m_sphere->get_color());
    
    m_sphere->draw();

    sp.unuse();
  }

  // WORLD COORDINATE AXES
  {
    auto& sp = m_sp_color_only;
    sp.use();
    sp.set_uniform("u_mvp", mvp);

    m_world_coord_axes->draw();

    sp.unuse();
  }
}