
#ifndef TRIANGLES_H
#define TRIANGLES_H

#include <vector>
#include <qvector3d.h>
#include "Common_defs.h"


class Triangles : protected OpenGLFunctionsBase
{
public:
  // IMPORTANT: we assume that the triangles are on the sphere!
  // this means that all vertex-normals are actually the normal vector on the
  // sphere at the point of projection of the current vertex on the sphere.
  Triangles(std::vector<QVector3D>& vertices);

  int get_num_triangles() const;

  void draw();
  
private:
  GLuint  m_vao, m_vbo;
  int     m_num_vertices;
};


#endif
