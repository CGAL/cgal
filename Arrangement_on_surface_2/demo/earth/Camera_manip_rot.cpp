
#include "Camera_manip_rot.h"


Camera_manip_rot::Camera_manip_rot(Camera& camera) : 
  Camera_manip(camera)
{
}

void Camera_manip_rot::mouse_press_event(QMouseEvent* e)
{
  set_mouse_button_pressed_flag(e, true);
  m_mouse_press_pos = m_last_mouse_pos = QVector2D(e->position());
}
void Camera_manip_rot::mouse_move_event(QMouseEvent* e)
{
  if (m_left_mouse_button_down)
  {
    auto current_mouse_pos = QVector2D(e->position());
    const auto diff = current_mouse_pos - m_last_mouse_pos;

    const float rotation_scale_factor = 0.1f;
    m_theta += rotation_scale_factor * diff.x();
    m_phi += rotation_scale_factor * diff.y();
    m_camera.rotate_from_init_config(-m_theta, -m_phi);
  }
}
void Camera_manip_rot::mouse_release_event(QMouseEvent* e)
{
}

