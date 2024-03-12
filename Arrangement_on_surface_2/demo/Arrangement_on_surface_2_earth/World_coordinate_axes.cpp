// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include <vector>

#include <qvector3d.h>

#include "World_coordinate_axes.h"

//! \brief
World_coord_axes::World_coord_axes(float length) {
  initializeOpenGLFunctions();

  const auto a = length;
  const float c = 0.0f;
  std::vector<QVector3D> vertex_data {
    QVector3D(0, 0, 0), QVector3D(1, 0, 0),
    QVector3D(a, 0, 0), QVector3D(1, 0, 0),

    QVector3D(0, 0, 0), QVector3D(0, 1, 0),
    QVector3D(0, a, 0), QVector3D(0, 1, 0),

    QVector3D(0, 0, 0), QVector3D(c, c, 1),
    QVector3D(0, 0, a), QVector3D(c, c, 1)
  };

  // DEFINE OPENGL BUFFERS
  glGenVertexArrays(1, &m_vao);
  glBindVertexArray(m_vao);

  // Vertex Buffer
  glGenBuffers(1, &m_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
  auto vertex_buffer_size = sizeof(QVector3D) * vertex_data.size();
  auto vertex_buffer_data = reinterpret_cast<const void*>(vertex_data.data());
  glBufferData(GL_ARRAY_BUFFER, vertex_buffer_size, vertex_buffer_data,
               GL_STATIC_DRAW);

  // Position Vertex-Attribute
  GLint position_attrib_index = 0;
  const void* position_offset = 0;
  GLsizei stride = 6 * sizeof(float);
  glVertexAttribPointer(position_attrib_index, 3, GL_FLOAT, GL_FALSE, stride,
                        position_offset);
  glEnableVertexAttribArray(position_attrib_index);

  // Color Vertex-Attribute
  GLint color_attrib_index = 1;
  auto* color_offset = reinterpret_cast<const void*>(3 * sizeof(float));
  glVertexAttribPointer(color_attrib_index, 3, GL_FLOAT, GL_FALSE, stride,
                        color_offset);
  glEnableVertexAttribArray(color_attrib_index);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

//! \brief
void World_coord_axes::draw() {
  glBindVertexArray(m_vao);
  const GLsizei count = 2 * 3; // = 2 * number of lines
  glDrawArrays(GL_LINES, 0, count);
  glBindVertexArray(0);
}
