
#include "Camera_manip_zoom.h"


Camera_manip_zoom::Camera_manip_zoom(Camera& camera) : 
  Camera_manip(camera)
{
}

void Camera_manip_zoom::mouse_move_event(QMouseEvent* e)
{
  if (m_middle_mouse_button_down)
  {
    auto current_mouse_pos = QVector2D(e->position());
    const auto diff = current_mouse_pos - m_last_mouse_pos;

    const float zoom_scale_factor = 0.01f;
    const auto distance = zoom_scale_factor * diff.y();
    m_camera.move_forward(distance);
  }
}

