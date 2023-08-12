// Copyright(c) 2012, 2020  Tel - Aviv University(Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include "Camera.h"


Camera::Camera() :
  m_ux(1, 0, 0),
  m_uy(0, 1, 0),
  m_uz(0, 0, 1)
{
}

void Camera::perspective(float fov, float aspect, float z_near, float z_far)
{
  m_z_near = z_near;
  m_z_far = z_far;

  m_projection.setToIdentity();
  m_projection.perspective(fov, aspect, z_near, z_far);
}


QMatrix4x4 Camera::get_view_matrix() const
{
  QMatrix4x4 view;
  const QVector3D center = m_pos - m_uz;
  view.lookAt(m_pos, center, m_uy);
  return view;
}


void Camera::rotate_from_init_config(float theta, float phi)
{
  // TO-DO: use the following logic to eliminate the QT-deprecation warnings!
  //QMatrix4x4 rot;
  //rot.rotate(theta, m_ux);
  //auto pos = m_pos.toVector4D();  pos.setW(1);
  //auto uy = m_uy.toVector4D();    uy.setW(0);
  //auto uz = m_uz.toVector4D();    uz.setW(0);

  //pos = pos * rot;
  //uy = uy * rot;
  //uz = uz * rot;

  //m_pos = pos.toVector3D();
  //m_uy = uy.toVector3D();
  //m_uz = uz.toVector3D();

  QMatrix4x4 r1;
  QVector3D ey(0, 1, 0);
  r1.rotate(theta, ey);

  // rx = rotated x axis
  auto rx = r1 * QVector3D(1,0,0);
  QMatrix4x4 r2;
  r2.rotate(phi, rx);

  // total rotation:
  auto r = r2 * r1;

  const auto dist_cam_to_origin = m_pos.length();
  m_pos = r * QVector3D(0, 0, dist_cam_to_origin);
  m_ux = r * QVector3D(1, 0, 0); // should be the same as rx (sanity check?)
  m_uy = r * QVector3D(0, 1, 0);
  m_uz = r * QVector3D(0, 0, 1);
}

void Camera::rotate(QMatrix4x4 rot)
{
  m_pos = rot * m_pos;
  m_ux = rot * m_ux;
  m_uy = rot * m_uy;
  m_uz = rot * m_uz;
}

void Camera::save_config()
{
  m_saved_pos = m_pos;
  m_saved_ux  = m_ux;
  m_saved_uy  = m_uy;
  m_saved_uz  = m_uz;
}

void Camera::rotate_from_saved_config(QMatrix4x4 rot)
{
  m_pos = rot * m_saved_pos;
  m_ux =  rot * m_saved_ux;
  m_uy =  rot * m_saved_uy;
  m_uz =  rot * m_saved_uz;
}

void Camera::move_forward(float distance)
{
  // recall that in OpenGL camera model, camera's z-axis points always
  // out of the screen (towards the user).
  m_pos -= distance * m_uz;
}
