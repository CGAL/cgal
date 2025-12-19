// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef GUI_COUNTRY_PICK_HANDLER_H
#define GUI_COUNTRY_PICK_HANDLER_H

#include <qevent.h>
#include <qvector2d.h>

#include "GUI_event_handler.h"
#include "Main_widget.h"

class GUI_country_pick_handler : public GUI_event_handler {
public:
  GUI_country_pick_handler(Main_widget& main_widget);

protected:
  virtual void mouse_press_event(QMouseEvent* e) override;
  virtual void resize(int w, int h) override;

  Main_widget& m_main_widget;
  Camera& m_camera;
  int m_vp_width;
  int m_vp_height;
};


#endif
