
#ifndef CAMERA_MANIP_ROT_H
#define CAMERA_MANIP_ROT_H

#include <qevent.h>
#include <qvector2d.h>

#include "Camera.h"


class Camera_manip_rot
{
public:
  Camera_manip_rot(Camera& camera);

  void mousePressEvent(QMouseEvent* e);
  void mouseMoveEvent(QMouseEvent* e);
  void mouseReleaseEvent(QMouseEvent* e);

protected:
  void set_mouse_button_pressed_flag(QMouseEvent* e, bool flag);

private:
  Camera& m_camera;
  float m_theta = 0, m_phi = 0;

  bool m_left_mouse_button_down = false;
  QVector2D m_last_mouse_pos;
  QVector2D m_mouse_press_pos;
};


#endif
