
#ifndef WORLD_COORD_AXES_H
#define WORLD_COORD_AXES_H

#include "Common_defs.h"
#include <qvector4d.h>


class World_coord_axes : protected OpenGLFunctionsBase
{
public:
  World_coord_axes(float length);

  void draw();

private:
  GLuint    m_vao, m_vbo;
};


#endif
