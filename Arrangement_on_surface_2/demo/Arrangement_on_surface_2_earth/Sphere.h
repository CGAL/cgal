// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef SPHERE_H
#define SPHERE_H

#include "Common_defs.h"
#include <qvector4d.h>

class Sphere : protected OpenGLFunctionsBase {
public:
  Sphere(std::size_t num_slices, std::size_t num_stacks, float r);

  void set_color(float r, float g, float b, float a);
  const QVector4D& get_color() { return m_color; }

  void draw();

private:
  GLuint m_vao;
  GLuint m_vbo;
  GLuint m_ibo;
  GLuint m_num_indices;
  QVector4D m_color;
};


#endif
