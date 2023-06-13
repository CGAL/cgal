
#ifndef SPHERE_H
#define SPHERE_H

#include "Common_defs.h"
#include <qvector4d.h>


class Sphere : protected OpenGLFunctionsBase
{
public:
  Sphere(int num_slices, int num_stacks, float r);

  void set_color(float r, float g, float b, float a);
  const QVector4D& get_color() { return m_color; }

  void draw();

private:
  GLuint    m_vao, m_vbo, m_ibo, m_num_indices;
  QVector4D m_color;
};


#endif
