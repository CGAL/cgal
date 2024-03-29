// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef CAMERA_H
#define CAMERA_H

#include <qvector3d.h>
#include <qmatrix4x4.h>

class Camera {
public:
  Camera();

  void set_pos(const QVector3D& pos) { m_pos = pos; }
  void set_pos(float x, float y, float z) { m_pos = QVector3D(x,y,z); }
  const QVector3D& get_pos() const { return m_pos; }

  void perspective(qreal fov, qreal aspect_ratio, qreal z_near, qreal z_far);

  qreal get_z_near() const { return m_z_near; }
  QMatrix4x4 get_view_matrix() const;
  QMatrix4x4 get_projection_matrix() const { return m_projection; }

  // theta: angle around y-axis
  // phi: angle from the xz-plane (= rotated x-axis after the above rotation)
  void rotate_from_init_config(float theta, float phi);
  void rotate(QMatrix4x4 rot);

  // save config & rotate from saved config (move to separate class?)
  void save_config();
  void rotate_from_saved_config(QMatrix4x4 rot);

  // move the camera forward around its own z-axis
  void move_forward(float distance);

private:
  QVector3D m_pos;
  QVector3D m_ux;
  QVector3D m_uy;
  QVector3D m_uz;

  QVector3D m_saved_pos;
  QVector3D m_saved_ux;
  QVector3D m_saved_uy;
  QVector3D m_saved_uz;

  qreal m_z_near, m_z_far;

  QMatrix4x4 m_projection;
};

#endif
