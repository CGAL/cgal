// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef VERTICES_H
#define VERTICES_H

#include <vector>
#include <qvector3d.h>

#include "Common_defs.h"

class Vertices : protected OpenGLFunctionsBase {
public:
  Vertices(const std::vector<QVector3D>& vertices);

  void draw();

private:
  GLuint m_vao;
  GLuint m_vbo;
  GLsizei m_num_indices;
  std::vector<GLuint> m_offsets;
};

#endif
