
#include "mainwidget.h"

#include <cmath>
#include <iostream>
#include <string>

#include <QMouseEvent>


MainWidget::~MainWidget()
{
  // Make sure the context is current when deleting the texture and the buffers.
  makeCurrent();
  doneCurrent();
}


float theta = 0, phi = 0;

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
void MainWidget::mousePressEvent(QMouseEvent* e)
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
    //const float theta_around_x = rotation_scale_factor * diff.y();
    //const float theta_around_y = rotation_scale_factor * diff.x();
    //m_camera.rotate(theta_around_x, theta_around_y);
    theta += rotation_scale_factor * diff.x();
    phi += rotation_scale_factor * diff.y();
    m_camera.rotate(-theta, -phi);
  }
  else if(m_middle_mouse_button_down)
  {
    const float zoom_scale_factor = 0.01f;
    const auto distance = zoom_scale_factor * diff.y();
    m_camera.move_forward(distance);
  }

  m_last_mouse_pos = current_mouse_pos;
}
void MainWidget::mouseReleaseEvent(QMouseEvent* e)
{
  set_mouse_button_pressed_flag(e, false);
}
void MainWidget::timerEvent(QTimerEvent*)
{
  update();
}


void MainWidget::initializeGL()
{
  initializeOpenGLFunctions();

  init_camera();
  init_geometry();
  init_shader_programs();

  {
    // has to be defined after camera has been defined:
    // because we want to compute the error based on camera parameters!
    m_geodesic_arcs = std::make_unique<Geodesic_arcs>();
  }

  glClearColor(0, 0, 0, 1);
  glEnable(GL_DEPTH_TEST);  // Enable depth buffer
  //glEnable(GL_CULL_FACE); // Enable back face culling

  // Use QBasicTimer because its faster than QTimer
  m_timer.start(12, this);
}



void MainWidget::init_camera()
{
  m_camera.set_pos(0, 0, 3);
  //m_camera.rotate_around_x(-90);
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
  init_sp_per_vertex_color();
  init_sp_arc();
}
void MainWidget::init_sp_smooth()
{
  const char* vs = "shaders/smooth_vs.glsl";
  const char* fs = "shaders/smooth_fs.glsl";
  m_sp_smooth.init(vs, "", fs);
}
void MainWidget::init_sp_per_vertex_color()
{
  const char* vs = "shaders/per_vertex_color_vs.glsl";
  const char* fs = "shaders/per_vertex_color_fs.glsl";
  m_sp_per_vertex_color.init(vs, "", fs);
}
void MainWidget::init_sp_arc()
{
  const char* vs = "shaders/arc_vs.glsl";
  const char* fs = "shaders/arc_fs.glsl";
  m_sp_arc.init(vs, "", fs);
}


std::ostream& operator << (std::ostream& os, const QVector4D& v)
{
  os << v.x() << ", " << v.y() << ", " << v.z() << ", " << v.w();
  return os;
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
    glEnable(GL_DEPTH_TEST);

    auto& sp = m_sp_smooth;
    sp.use();
    sp.set_uniform("u_mvp", mvp);
    sp.set_uniform("u_color", m_sphere->get_color());
    
    m_sphere->draw();

    sp.unuse();
  }

  // WORLD COORDINATE AXES &  GEODESIC ARCS
  {
    auto& sp = m_sp_per_vertex_color;
    sp.use();
    sp.set_uniform("u_mvp", mvp);

    m_world_coord_axes->draw();

    sp.unuse();
  }

  // GEODESIC ARCS
  {
    glDisable(GL_DEPTH_TEST);

    auto& sp = m_sp_arc;
    sp.use();
    sp.set_uniform("u_mvp", mvp);

    // compute the cutting plane
    auto c = m_camera.get_pos();
    const auto d = c.length();
    const auto r = 1.0f;
    const auto sin_alpha = r / d;
    const auto n = (c / d); // plane unit normal vector
    const auto cos_beta = sin_alpha;
    const auto p = (r * cos_beta) * n;
    QVector4D plane(n.x(), n.y(), n.z(), -QVector3D::dotProduct(p, n));
    const QVector4D arc_color(0, 0.5, 1, 1);
    glLineWidth(5);
    sp.set_uniform("u_color", arc_color);
    sp.set_uniform("u_plane", plane);
    m_geodesic_arcs->draw();

    sp.unuse();
  }
}