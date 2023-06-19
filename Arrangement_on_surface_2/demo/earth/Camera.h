
#ifndef CAMERA_H
#define CAMERA_H


#include <qvector3d.h>
#include <qmatrix4x4.h>


class Camera
{
public:

  Camera();

  void set_pos(const QVector3D& pos) { m_pos = pos; }
  void set_pos(float x, float y, float z) { m_pos = QVector3D(x,y,z); }
  const QVector3D& get_pos() const { return m_pos; }

  void perspective(float fov, float aspect_ratio, float z_near, float z_far);

  float get_z_near() const { return m_z_near; }
  QMatrix4x4 get_view_matrix() const;
  QMatrix4x4 get_projection_matrix() const { return m_projection; }
  
  // theta: angle around y-axis
  // phi: angle from the xz-plane (= rotated x-axis after the above rotation)
  void rotate(float theta, float phi);
  void rotate(QMatrix4x4 rot);

  // move the camera forward around its own z-axis
  void move_forward(float distance);

private:
  QVector3D m_pos;
  QVector3D m_ux;
  QVector3D m_uy;
  QVector3D m_uz;

  float m_z_near, m_z_far;

  QMatrix4x4 m_projection;
};


#endif
