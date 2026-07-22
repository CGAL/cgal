// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef LINE_STRIPS_H
#define LINE_STRIPS_H

#include <vector>
#include <qvector3d.h>

#include "Common_defs.h"

class Line_strips : protected OpenGLFunctionsBase {
public:
  Line_strips(std::vector<QVector3D>& line_strip_points);
  Line_strips(std::vector<std::vector<QVector3D>>& arcs);

  std::size_t get_num_line_strips() const;
  void draw(int line_strip_index);

  void draw();

private:
  GLuint m_vao;
  GLuint m_vbo;
  std::vector<GLuint> m_offsets;
};


#endif
