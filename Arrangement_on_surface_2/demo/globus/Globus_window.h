// Copyright (c) 2022 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel <efifogel@gmail.com>

#ifndef GLOBUS_WINDOW_H
#define GLOBUS_WINDOW_H

#include <utility>
#include <vector>

#include <Qt>

#include <CGAL/Object.h>
#include <CGAL/Qt/DemosMainWindow.h>

namespace Ui { class Globus_window; }

namespace CGAL {
class Object;
namespace Qt {
class GraphicsViewNavigation;
}
}

class QActionGroup;

class Globus_window : public CGAL::Qt::DemosMainWindow {
  Q_OBJECT

public:
  Globus_window(QWidget* parent = nullptr);
  ~Globus_window();

public Q_SLOTS:
  void on_action_quit_triggered();
  void on_action_open_triggered();
  void on_action_save_as_triggered();
  void on_action_preferences_triggered();
  void on_action_zoom_in_triggered();
  void on_action_zoom_out_triggered();
  void on_action_zoom_reset_triggered();
  void on_action_drag_toggled(bool);

Q_SIGNALS:

protected:
  void setup_ui();

private:
  Ui::Globus_window* ui;
  QActionGroup* mode_group;
};

#endif
