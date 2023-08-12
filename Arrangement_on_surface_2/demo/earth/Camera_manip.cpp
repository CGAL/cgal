// Copyright(c) 2012, 2020  Tel - Aviv University(Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include "Camera_manip.h"


Camera_manip::Camera_manip(Camera& camera) : 
  m_camera(camera)
{
}


void Camera_manip::set_mouse_button_pressed_flag(QMouseEvent* e, bool flag)
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
void Camera_manip::mousePressEvent(QMouseEvent* e)
{
  set_mouse_button_pressed_flag(e, true);
  m_mouse_press_pos = m_last_mouse_pos = QVector2D(e->position());
  
  // call the function overridden by the derived class
  mouse_press_event(e);
}
void Camera_manip::mouseMoveEvent(QMouseEvent* e)
{
  m_current_mouse_pos = QVector2D(e->position());
  m_diff = m_current_mouse_pos - m_last_mouse_pos;

  // call the function overridden by the derived class
  mouse_move_event(e);

  m_last_mouse_pos = m_current_mouse_pos;
}
void Camera_manip::mouseReleaseEvent(QMouseEvent* e)
{
  set_mouse_button_pressed_flag(e, false);

  // call the function overridden by the derived class
  mouse_release_event(e);
}
void Camera_manip::resizeGL(int w, int h)
{
  resize(w, h);
}

