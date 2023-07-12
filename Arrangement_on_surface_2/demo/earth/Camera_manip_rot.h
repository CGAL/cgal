
#ifndef CAMERA_MANIP_ROT_H
#define CAMERA_MANIP_ROT_H

#include <qevent.h>
#include <qvector2d.h>

#include "Camera_manip.h"


class Camera_manip_rot : public Camera_manip
{
public:
  Camera_manip_rot(Camera& camera);

protected:
  void mouse_press_event(QMouseEvent* e) override;
  void mouse_move_event(QMouseEvent* e) override;
  void mouse_release_event(QMouseEvent* e) override;

private:
  float m_theta = 0, m_phi = 0;

  //bool m_left_mouse_button_down = false;
  //QVector2D m_last_mouse_pos;
  //QVector2D m_mouse_press_pos;
};


#endif
