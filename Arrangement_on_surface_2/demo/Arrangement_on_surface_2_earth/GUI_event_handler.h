// Copyright (c) 2023, 2024  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef GUI_EVENT_HANDLER_H
#define GUI_EVENT_HANDLER_H

#include <qevent.h>
#include <qvector2d.h>

class GUI_event_handler {
public:
  virtual ~GUI_event_handler() {};

  void mousePressEvent(QMouseEvent* e);
  void mouseMoveEvent(QMouseEvent* e);
  void mouseReleaseEvent(QMouseEvent* e);
  void resizeGL(int w, int h);

protected:
  void set_mouse_button_pressed_flag(QMouseEvent* e, bool flag);

  virtual void mouse_press_event(QMouseEvent* /* e */) {}
  virtual void mouse_move_event(QMouseEvent* /* e */) {}
  virtual void mouse_release_event(QMouseEvent* /* e */) {}
  virtual void resize(int /* w */, int /* h */) {}

  bool m_left_mouse_button_down = false;
  bool m_middle_mouse_button_down = false;
  QVector2D m_current_mouse_pos;
  QVector2D m_last_mouse_pos;
  QVector2D m_mouse_press_pos;
  QVector2D m_diff;
};


#endif
