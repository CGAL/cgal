// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include "Camera_manip_rot.h"

Camera_manip_rot::Camera_manip_rot(Camera& camera) : Camera_manip(camera) {}

void Camera_manip_rot::mouse_move_event(QMouseEvent* /* e */) {
  if (m_left_mouse_button_down) {
    const float rotation_scale_factor = 0.1f;
    m_theta += rotation_scale_factor * m_diff.x();
    m_phi += rotation_scale_factor * m_diff.y();
    m_camera.rotate_from_init_config(-m_theta, -m_phi);
  }
}
