
#ifndef GEODESIC_ARCS_H
#define GEODESIC_ARCS_H

#include "Common_defs.h"


class Geodesic_arcs : protected OpenGLFunctionsBase
{
public:
  Geodesic_arcs();

  void draw();

private:
  GLuint    m_vao, m_vbo, m_num_arc_points;
};


#endif
