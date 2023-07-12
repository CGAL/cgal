
#include "Camera_manip_rot.h"


Camera_manip_rot::Camera_manip_rot(Camera& camera) : 
  m_camera(camera)
{
}


void Camera_manip_rot::set_mouse_button_pressed_flag(QMouseEvent* e, bool flag)
{
  switch (e->button())
  {
  case Qt::LeftButton:
    m_left_mouse_button_down = flag;
    break;
  }
}
void Camera_manip_rot::mousePressEvent(QMouseEvent* e)
{
  set_mouse_button_pressed_flag(e, true);
  m_mouse_press_pos = m_last_mouse_pos = QVector2D(e->position());
}
void Camera_manip_rot::mouseMoveEvent(QMouseEvent* e)
{
  auto current_mouse_pos = QVector2D(e->position());
  const auto diff = current_mouse_pos - m_last_mouse_pos;

  if (m_left_mouse_button_down)
  {
    const float rotation_scale_factor = 0.1f;
    m_theta += rotation_scale_factor * diff.x();
    m_phi += rotation_scale_factor * diff.y();
    m_camera.rotate_from_init_config(-m_theta, -m_phi);
  }

  m_last_mouse_pos = current_mouse_pos;
}
void Camera_manip_rot::mouseReleaseEvent(QMouseEvent* e)
{
  set_mouse_button_pressed_flag(e, false);
}

