
#ifndef LINE_STRIPS_H
#define LINE_STRIPS_H

#include <vector>
#include "Common_defs.h"


class Line_strips : protected OpenGLFunctionsBase
{
public:
  Line_strips(Line_strip_approx& lsa);

  void draw();

private:
  GLuint                m_vao, m_vbo;
  std::vector<GLuint>   m_offsets;
};


#endif
