
#ifndef GEODESIC_ARCS_H
#define GEODESIC_ARCS_H

#include "Common_defs.h"
#include <vector>

class Geodesic_arcs : protected OpenGLFunctionsBase
{
public:
  Geodesic_arcs();

  void draw();

private:
  GLuint    m_vao, m_vbo, m_num_arc_points;
  std::vector<GLuint>   m_arc_offsets;
};


#endif
