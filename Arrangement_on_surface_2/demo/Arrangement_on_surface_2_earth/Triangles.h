// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef TRIANGLES_H
#define TRIANGLES_H

#include <vector>
#include <qvector3d.h>
#include "Common_defs.h"

class Triangles : protected OpenGLFunctionsBase {
public:
  // IMPORTANT: we assume that the triangles are on the sphere!
  // this means that all vertex-normals are actually the normal vector on the
  // sphere at the point of projection of the current vertex on the sphere.
  Triangles(std::vector<QVector3D>& vertices);

  int get_num_triangles() const;
  void set_color(const QVector4D& rgba);
  const QVector4D& get_color() const;

  void draw();

private:
  GLuint m_vao;
  GLuint m_vbo;
  GLsizei m_num_vertices;
  QVector4D m_color;
};


#endif
