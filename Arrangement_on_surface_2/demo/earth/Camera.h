
#ifndef CAMERA_H
#define CAMERA_H


#include <qvector3d.h>
#include <qmatrix4x4.h>


class Camera
{
public:

  Camera() :
    m_ux(1, 0, 0),
    m_uy(0, 1, 0),
    m_uz(0, 0, 1)
  {
  }

  void set_pos(const QVector3D& pos) { m_pos = pos; }
  void set_pos(float x, float y, float z) { m_pos = QVector3D(x,y,z); }

  QMatrix4x4 get_view_matrix() const
  {
    QMatrix4x4 view;
    const QVector3D center = m_pos - m_uz;
    view.lookAt(m_pos, center, m_uy);
    return view;
  }

  // rotate the camera around its own axes
  void rotate(float theta_around_x, float theta_around_y)
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

private:
  QVector3D m_pos;
  QVector3D m_ux;
  QVector3D m_uy;
  QVector3D m_uz;
};


#endif
