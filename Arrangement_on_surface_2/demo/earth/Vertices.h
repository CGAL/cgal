
#ifndef VERTICES_H
#define VERTICES_H

#include <vector>
#include <qvector3d.h>
#include "Common_defs.h"


class Vertices : protected OpenGLFunctionsBase
{
public:
  Vertices(const std::vector<QVector3D>& vertices);

  void draw();

private:
  GLuint                m_vao, m_vbo, m_num_indices;
  std::vector<GLuint>   m_offsets;
};


#endif
