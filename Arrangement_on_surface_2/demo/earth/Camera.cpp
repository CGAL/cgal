
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


void Camera::rotate(float theta_around_x, float theta_around_y)
{
  // rotate the camera around its x-axis
  QMatrix4x4 rot;
  rot.rotate(theta_around_x, m_ux);
  m_pos = m_pos * rot;
  m_uy = m_uy * rot;
  m_uz = m_uz * rot;

  // rotate the camera around its y-axis
  rot.setToIdentity();
  rot.rotate(theta_around_y, m_uy);
  m_pos = m_pos * rot;
  m_ux = m_ux * rot;
  m_uz = m_uz * rot;
}

void Camera::move_forward(float distance)
{
  // recall that in OpenGL camera model, camera's z-axis points always
  // out of the screen (towards the user).
  m_pos -= distance * m_uz;
}
