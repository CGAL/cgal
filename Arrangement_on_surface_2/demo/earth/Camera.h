
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

  QMatrix4x4 get_view_matrix() const;
  QMatrix4x4 get_projection_matrix() const { return m_projection; }

  // rotate the camera around its own axes
  void rotate_around_x(float theta);
  void rotate_around_y(float theta);
  void rotate(float theta_around_x, float theta_around_y);

  // move the camera forward around its own z-axis
  void move_forward(float distance);

private:
  QVector3D m_pos;
  QVector3D m_ux;
  QVector3D m_uy;
  QVector3D m_uz;

  QMatrix4x4 m_projection;
};


#endif
