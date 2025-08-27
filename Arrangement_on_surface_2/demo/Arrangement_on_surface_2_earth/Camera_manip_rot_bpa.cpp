// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include "Camera_manip_rot_bpa.h"

//! \brief
Camera_manip_rot_bpa::Camera_manip_rot_bpa(Camera& camera) :
  Camera_manip(camera)
{}

//! \brief
void Camera_manip_rot_bpa::mouse_press_event(QMouseEvent* /* e */) {
  // for the backprojected diff-vector method:
  if (m_left_mouse_button_down) m_camera.save_config();
}

//! \brief
void Camera_manip_rot_bpa::mouse_move_event(QMouseEvent* /* e */) {
  const float rotation_scale_factor = 0.1f;

  // ROTATION AROUND AN AXIS ORTHOGONAL TO THE BACKPROJECTED DIF-VECTOR
  //QVector3D p0(m_last_mouse_pos.x(), m_vp_height - m_last_mouse_pos.y(), 0);
  QVector3D p0(m_mouse_press_pos.x(), m_vp_height - m_mouse_press_pos.y(), 0);
  QVector3D p1(m_current_mouse_pos.x(), m_vp_height - m_current_mouse_pos.y(), 0);
  auto dp = p1 - p0; // difference vector in OpenGL window coords.
  QVector3D rdp(-dp.y(), dp.x(), 0); // rotate diff-vector CCW by 90-deg
  QVector3D rp = p0 + rdp; // r1 rotated CCW by 90 deg

  QMatrix4x4 model; // this is different from Sphere's model matrix!!!
  auto proj = m_camera.get_projection_matrix();
  auto view = m_camera.get_view_matrix();
  auto model_view = view * model;
  QRect viewport(0, 0, m_vp_width, m_vp_height);
  auto wp0 = p0.unproject(model_view, proj, viewport);
  auto wrp = rp.unproject(model_view, proj, viewport);

  // rotation axis & angle
  auto rot_axis = wrp - wp0;
  rot_axis.normalize();
  const auto rot_angle = rotation_scale_factor * dp.length();

  QMatrix4x4 rot_matrix;
  rot_matrix.rotate(-rot_angle, rot_axis);

  m_camera.rotate_from_saved_config(rot_matrix);
}

//! \brief
void Camera_manip_rot_bpa::mouse_release_event(QMouseEvent* /* e */) {}

//! \brief
void Camera_manip_rot_bpa::resize(int w, int h) {
  m_vp_width = w;
  m_vp_height = h;
}
