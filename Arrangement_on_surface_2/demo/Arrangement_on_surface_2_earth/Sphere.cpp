// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include <cmath>
#include <vector>

#include <qvector3d.h>

#include "Sphere.h"

//! \brief
Sphere::Sphere(std::size_t num_slices, std::size_t num_stacks, float r) {
  initializeOpenGLFunctions();

  num_stacks = std::max<std::size_t>(2, num_stacks);
  std::vector<QVector3D> vertices, normals;

  // NORTH POLE
  vertices.push_back(QVector3D(0, 0, r));
  normals.push_back(QVector3D(0, 0, 1));

  // SOUTH POLE
  vertices.push_back(QVector3D(0, 0, -r));
  normals.push_back(QVector3D(0, 0, -1));
  auto starting_index_of_middle_vertices = vertices.size();

  for (std::size_t j = 1; j < num_stacks; ++j) {
    // Calculate the latitude (vertical angle) for the current stack
    float lat = M_PI * j / num_stacks;
    float rxy = r * std::sin(lat);
    float z = r * std::cos(lat);

    for (std::size_t i = 0; i < num_slices; ++i) {
      // Calculate the longitude (horizontal angle) for the current slice
      float lon = 2 * M_PI * i / num_slices;

      // Convert spherical coordinates to Cartesian coordinates
      float x = rxy * std::cos(lon);
      float y = rxy * std::sin(lon);

      auto p = QVector3D(x, y, z);
      auto n = p / p.length();
      vertices.push_back(p);
      normals.push_back(n);
    }
  }

  // strided vertex-data
  std::vector<QVector3D> vertex_data;
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    vertex_data.push_back(vertices[i]);
    vertex_data.push_back(normals[i]);
  }

  // add the indices for all triangles
  std::vector<GLuint> indices;

  // NORTH CAP
  const GLuint north_vertex_index = 0;
  const auto north_cap_vertex_index_start = starting_index_of_middle_vertices;
  for (std::size_t i = 0; i < num_slices; ++i) {
    indices.push_back(north_vertex_index);
    indices.push_back(static_cast<GLuint>(north_cap_vertex_index_start + i));
    indices.push_back(static_cast<GLuint>(north_cap_vertex_index_start + (i + 1) % num_slices));
  }

  // 0 = NORTH VERTEX
  // 1 = SOUTH VERTEX
  // [2, 2 + (numSlices-1)] = bottom vertices of the stack #1
  // [2+numSlices, 2 + (2*numSlices - 1)] = bottom vertices of the stack #2
  // ...
  // [2+(k-1)*numSlices, 2 + (k*numSlices -1)] = bottom vertices of the stack #k
  // ..
  // [2+(numStacks-1)*numSlices, 2+(numStacks*numSlices-1)] = bottom vertices of
  //                                                the last stack (# numStacks)

  // SOUTH CAP
  const GLuint south_vertex_index = 1;
  const std::size_t south_cap_index_start = starting_index_of_middle_vertices +
    (num_stacks - 2) * num_slices;
  for (std::size_t i = 0; i < num_slices; ++i) {
    const auto vi0 = south_vertex_index;
    const auto vi1 = static_cast<GLuint>(south_cap_index_start + i);
    const auto vi2 = static_cast<GLuint>(south_cap_index_start + (i + 1) % num_slices);
    indices.push_back(vi2);
    indices.push_back(vi1);
    indices.push_back(vi0);
  }

  // MIDDLE TRIANGLES
  for (std::size_t k = 0; k < num_stacks - 2; ++k) {
    const std::size_t stack_start_index =
      starting_index_of_middle_vertices + k * num_slices;
    const std::size_t next_stack_start_index = stack_start_index + num_slices;
    for (std::size_t i = 0; i < num_slices; ++i) {
      // check why the following code snippet does not work (winding order?)
      //std::size_t vi0 = stackStartIndex + i;
      //std::size_t vi1 = nextStackStartIndex + i;
      //std::size_t vi2 = nextStackStartIndex + (i + 1) % numSlices;
      //std::size_t vi3 = stackStartIndex + (i + 1) % numSlices;
      auto vi0 = static_cast<GLuint>(stack_start_index + i);
      auto vi1 = static_cast<GLuint>(stack_start_index + (i + 1) % num_slices);
      auto vi2 = static_cast<GLuint>(next_stack_start_index + i);
      auto vi3 = static_cast<GLuint>(next_stack_start_index + (i + 1) % num_slices);

      indices.push_back(vi0);
      indices.push_back(vi2);
      indices.push_back(vi1);
      //
      indices.push_back(vi2);
      indices.push_back(vi3);
      indices.push_back(vi1);
    }
  }
  m_num_indices = static_cast<GLuint>(indices.size());

  // DEFINE OPENGL BUFFERS
  glGenVertexArrays(1, &m_vao);
  glBindVertexArray(m_vao);

  // Index buffer
  glGenBuffers(1, &m_ibo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
  auto indices_size = sizeof(GLuint) * indices.size();
  auto indices_data = reinterpret_cast<const void*>(indices.data());
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices_size, indices_data,
               GL_STATIC_DRAW);

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

  // Note: calling this before glBindVertexArray(0) results in no output!
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

//! \brief
void Sphere::set_color(float r, float g, float b, float a)
{ m_color = QVector4D(r, g, b, a); }

//! \brief
void Sphere::draw() {
  // DRAW TRIANGLES
  glBindVertexArray(m_vao);

  //glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
  glDrawElements(GL_TRIANGLES, m_num_indices, GL_UNSIGNED_INT, 0);

  glBindVertexArray(0);
}
