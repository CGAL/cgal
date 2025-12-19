// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include "Triangles.h"

#include <iostream>

//! \brief
Triangles::Triangles(std::vector<QVector3D>& vertices) {
  initializeOpenGLFunctions();

  // computer the normals of each vertex
  std::vector<QVector3D> normals;
  for(const auto& p : vertices) {
    auto n = p;
    n.normalize();
    normals.push_back(n);
  }

  // std::size_t num_triangles = vertices.size() / 3;

  // strided vertex-data
  std::vector<QVector3D> vertex_data;
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    vertex_data.push_back(vertices[i]);
    vertex_data.push_back(normals[i]);
  }

  // DEFINE OPENGL BUFFERS
  glGenVertexArrays(1, &m_vao);
  glBindVertexArray(m_vao);
  m_num_vertices = static_cast<GLsizei>(vertices.size());

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

  // Normal Vertex-Attribute
  GLint normal_attrib_index = 1;
  auto* normal_offset = reinterpret_cast<const void*>(3 * sizeof(float));
  glVertexAttribPointer(normal_attrib_index, 3, GL_FLOAT, GL_FALSE, stride,
                        normal_offset);
  glEnableVertexAttribArray(normal_attrib_index);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

//! \brief
void Triangles::set_color(const QVector4D& rgba)
{ m_color = rgba; }

//! \brief
const QVector4D& Triangles::get_color() const { return m_color; }

//! \brief
void Triangles::draw() {
  // DRAW TRIANGLES
  glBindVertexArray(m_vao);
  glDrawArrays(GL_TRIANGLES, 0, m_num_vertices);
  glBindVertexArray(0);
}
