
#ifndef LINE_STRIPS_H
#define LINE_STRIPS_H

#include <vector>
#include <qvector3d.h>
#include "Common_defs.h"


class Line_strips : protected OpenGLFunctionsBase
{
public:
  Line_strips(std::vector<std::vector<QVector3D>>& arcs);

  void draw();

private:
  GLuint                m_vao, m_vbo;
  std::vector<GLuint>   m_offsets;
};


#endif
