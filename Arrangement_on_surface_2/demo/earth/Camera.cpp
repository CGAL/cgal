
#include "Camera.h"


Camera::Camera() :
  m_ux(1, 0, 0),
  m_uy(0, 1, 0),
  m_uz(0, 0, 1)
{
}

void Camera::perspective(float fov, float aspect, float z_near, float z_far)
{
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


void Camera::rotate_around_x(float theta)
{
  QMatrix4x4 rot;
  rot.rotate(theta, m_ux);
  auto pos = m_pos.toVector4D();  pos.setW(1);
  auto uy = m_uy.toVector4D();    uy.setW(0);
  auto uz = m_uz.toVector4D();    uz.setW(0);

  pos = pos * rot;
  uy = uy * rot;
  uz = uz * rot;

  m_pos = pos.toVector3D();
  m_uy = uy.toVector3D();
  m_uz = uz.toVector3D();
}
void Camera::rotate_around_y(float theta)
{
  QMatrix4x4 rot;
  rot.rotate(theta, m_uy);
  auto pos = m_pos.toVector4D();  pos.setW(1);
  auto ux = m_ux.toVector4D();    ux.setW(0);
  auto uz = m_uz.toVector4D();    uz.setW(0);

  pos = pos * rot;
  ux = ux * rot;
  uz = uz * rot;

  m_pos = pos.toVector3D();
  m_ux = ux.toVector3D();
  m_uz = uz.toVector3D();
}
void Camera::rotate(float theta_around_x, float theta_around_y)
{
  rotate_around_x(theta_around_x);
  rotate_around_y(theta_around_y);
}

void Camera::move_forward(float distance)
{
  // recall that in OpenGL camera model, camera's z-axis points always
  // out of the screen (towards the user).
  m_pos -= distance * m_uz;
}
