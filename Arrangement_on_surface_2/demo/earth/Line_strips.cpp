
#include "Line_strips.h"


Line_strips::Line_strips(std::vector<std::vector<QVector3D>>& arcs)
{
  initializeOpenGLFunctions();

  std::vector<QVector3D> vertex_data;
  m_offsets.push_back(0);
  for (const auto& arc : arcs)
  {
    for(const auto& p : arc)
      vertex_data.push_back(p);
    
    const auto end_of_current_arc_points = vertex_data.size();
    m_offsets.push_back(end_of_current_arc_points);
  }


  // DEFINE OPENGL BUFFERS
  glGenVertexArrays(1, &m_vao);
  glBindVertexArray(m_vao);


  // Vertex Buffer
  glGenBuffers(1, &m_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
  auto vertex_buffer_size = sizeof(QVector3D) * vertex_data.size();
  auto vertex_buffer_data = reinterpret_cast<const void*>(vertex_data.data());
  glBufferData(GL_ARRAY_BUFFER,
    vertex_buffer_size,
    vertex_buffer_data,
    GL_STATIC_DRAW);

  // Position Vertex-Attribute
  GLint position_attrib_index = 0;
  const void* position_offset = 0;
  GLsizei stride = 0;
  glVertexAttribPointer(position_attrib_index,
                        3,
                        GL_FLOAT, GL_FALSE,
                        stride,
                        position_offset);
  glEnableVertexAttribArray(position_attrib_index);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

void Line_strips::draw()
{
  glBindVertexArray(m_vao);
  {
    for (int i = 1; i < m_offsets.size(); i++)
    {
      const auto first = m_offsets[i - 1];
      const auto count = m_offsets[i] - first;
      glDrawArrays(GL_LINE_STRIP, first, count);
    }
  }
  glBindVertexArray(0);
}