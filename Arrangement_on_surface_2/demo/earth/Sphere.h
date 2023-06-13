
#ifndef SPHERE_H
#define SPHERE_H

#include "Common_defs.h"


class Sphere : protected OpenGLFunctionsBase
{
public:
  Sphere(int num_slices, int num_stacks, float r);

  void draw();

private:
  GLuint m_vao, m_vbo, m_ibo, m_num_indices;
};


#endif
