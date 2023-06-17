
#include "Main_widget.h"

#include <cmath>
#include <iostream>
#include <string>

#include <QMouseEvent>

#include "Geodesic_arcs.h"
#include "Tools.h"


Main_widget::~Main_widget()
{
  // Make sure the context is current when deleting the texture and the buffers.
  makeCurrent();
  doneCurrent();
}

void Main_widget::set_mouse_button_pressed_flag(QMouseEvent* e, bool flag)
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
void Main_widget::mousePressEvent(QMouseEvent* e)
{
  set_mouse_button_pressed_flag(e, true);
  m_last_mouse_pos = QVector2D(e->position());
}
void Main_widget::mouseMoveEvent(QMouseEvent* e)
{
  auto current_mouse_pos = QVector2D(e->position());
  const auto diff = current_mouse_pos - m_last_mouse_pos;

  if (m_left_mouse_button_down)
  {
    const float rotation_scale_factor = 0.1f;

    if(1)
    {
      // OUR CUSTOM AD-HOC CAMERA ROTATION
      m_theta += rotation_scale_factor * diff.x();
      m_phi += rotation_scale_factor * diff.y();
      m_camera.rotate(-m_theta, -m_phi);
    }
    else
    {
      // ROTATION AROUND AN AXIS ORTHOGONAL TO THE BACKPROJECTED DIF-VECTOR
      QVector3D p0(m_last_mouse_pos.x(), m_vp_height - m_last_mouse_pos.y(), 0);
      QVector3D p1(current_mouse_pos.x(), m_vp_height - current_mouse_pos.y(), 0);
      auto dp = p1 - p0; // difference vector in OpenGL window coords.
      QVector3D rdp(-dp.y(), dp.x(), 0); // rotate diff-vector CCW by 90-deg
      QVector3D rp = p0 + rdp; // r1 rotated CCW by 90 deg
     
      QMatrix4x4 model; // this is different from Sphere's model matrix!!!
      auto proj = m_camera.get_projection_matrix();
      auto view = m_camera.get_view_matrix();
      auto model_view = view * model;
      QRect viewport(0, 0, m_vp_width, m_vp_height);
      auto wp0 = p0.unproject(model_view, proj, viewport);
      auto wrp = rp.unproject(model_view, proj, viewport);

      // rotation axis & angle
      auto rot_axis = wrp - wp0;
      rot_axis.normalize();
      const auto rot_angle = rotation_scale_factor * dp.length();
    
      QMatrix4x4 rot_matrix;
      rot_matrix.rotate(-rot_angle, rot_axis);
      
      m_camera.rotate(rot_matrix);
    }
  }
  else if(m_middle_mouse_button_down)
  {
    const float zoom_scale_factor = 0.01f;
    const auto distance = zoom_scale_factor * diff.y();
    m_camera.move_forward(distance);
  }

  m_last_mouse_pos = current_mouse_pos;
}
void Main_widget::mouseReleaseEvent(QMouseEvent* e)
{
  set_mouse_button_pressed_flag(e, false);
}
void Main_widget::timerEvent(QTimerEvent*)
{
  update();
}


void Main_widget::initializeGL()
{
  initializeOpenGLFunctions();

  init_camera();
  init_geometry();
  init_shader_programs();

  {
    // has to be defined after camera has been defined:
    // because we want to compute the error based on camera parameters!
    Geodesic_arcs ga;
    const double error = 0.001; // calculate this from cam parameters!
    auto lsa = ga.get_approximate_arcs(error);
    m_geodesic_arcs = std::make_unique<Line_strips>(lsa);
  }

  glClearColor(0, 0, 0, 1);
  glEnable(GL_DEPTH_TEST);  // Enable depth buffer
  //glEnable(GL_CULL_FACE); // Enable back face culling

  // Use QBasicTimer because its faster than QTimer
  m_timer.start(12, this);
}



void Main_widget::init_camera()
{
  m_camera.set_pos(0, 0, 3);
  //m_camera.rotate_around_x(-90);
}
void Main_widget::init_geometry()
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
void Main_widget::init_shader_programs()
{
  Shader_program::set_shader_path("shaders/");
  m_sp_smooth.init_with_vs_fs("smooth");;
  m_sp_per_vertex_color.init_with_vs_fs("per_vertex_color");
  m_sp_arc.init_with_vs_fs("arc");
}

void Main_widget::resizeGL(int w, int h)
{
  m_vp_width = w;
  m_vp_height = h;

  // Reset projection
  qreal aspect = qreal(w) / qreal(h ? h : 1);
  const qreal z_near = 1.0, z_far = 100.0, fov = 45.0;
  m_camera.perspective(fov, aspect, z_near, z_far);
}
void Main_widget::paintGL()
{
  QMatrix4x4 model;
  model.rotate(-90, 1,0,0); // this makes z-axes point upwards!
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
    //glDisable(GL_DEPTH_TEST);

    auto& sp = m_sp_arc;
    sp.use();
    sp.set_uniform("u_mvp", mvp);

    // compute the cutting plane
    // remember that we are passing the local vertex positions of the sphere 
    // between the vertex and fragment shader stages, so we need to convert
    // the camera-pos in world coords to sphere's local coords!
    auto c =  model.inverted() * m_camera.get_pos();
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
