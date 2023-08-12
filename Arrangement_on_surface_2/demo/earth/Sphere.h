// Copyright(c) 2012, 2020  Tel - Aviv University(Israel).
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


class Sphere : protected OpenGLFunctionsBase
{
public:
  Sphere(int num_slices, int num_stacks, float r);

  void set_color(float r, float g, float b, float a);
  const QVector4D& get_color() { return m_color; }

  void draw();

private:
  GLuint    m_vao, m_vbo, m_ibo, m_num_indices;
  QVector4D m_color;
};


#endif
