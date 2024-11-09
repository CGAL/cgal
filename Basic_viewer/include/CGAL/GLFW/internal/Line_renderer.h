#ifndef CGAL_GLFW_INTERNAL_LINE_RENDERER_H
#define CGAL_GLFW_INTERNAL_LINE_RENDERER_H

#include <vector>

#include <glad/glad.h>

#include "utils.h"

class Line_renderer {
public:
  ~Line_renderer();

  void initialize_buffers();
  void load_buffers();
  void add_line(const vec3f& start, const vec3f& end);
  void add_line(const vec3f& start, const vec3f& end, const vec3f& color);
  void draw();

  inline bool are_buffers_loaded() const { return m_AreBuffersLoaded; }
  inline bool are_buffers_initialized() const { return m_AreBuffersInitialized; }

  inline void set_width(const float width) { m_Width = width; }

private: 
  std::vector<float> m_Vertices {};

  float m_Width { 1.0f };

  unsigned int m_VertexArray;
  unsigned int m_VertexBuffer;

  bool m_AreBuffersLoaded { false };
  bool m_AreBuffersInitialized { false };
};

/********************METHOD IMPLEMENTATIONS********************/

Line_renderer::~Line_renderer() 
{
  if (m_VertexArray != 0) glDeleteVertexArrays(1, &m_VertexArray);
  if (m_VertexBuffer != 0) glDeleteBuffers(1, &m_VertexBuffer);
}

void Line_renderer::initialize_buffers() 
{
  glGenVertexArrays(1, &m_VertexArray);
  glGenBuffers(1, &m_VertexBuffer);

  glBindVertexArray(m_VertexArray);

  glBindBuffer(GL_ARRAY_BUFFER, m_VertexBuffer);

  // vertex attribute  
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), nullptr);

  // color attribute 
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (const void*)12);

  glBindVertexArray(0);  
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  m_AreBuffersInitialized = true;
}

void Line_renderer::load_buffers() 
{
  assert(are_buffers_initialized());

  glBindVertexArray(m_VertexArray);

  glBindBuffer(GL_ARRAY_BUFFER, m_VertexBuffer);
  glBufferData(GL_ARRAY_BUFFER, m_Vertices.size() * sizeof(float), m_Vertices.data(), GL_STATIC_DRAW);

  glBindVertexArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  m_AreBuffersLoaded = true;
}

void Line_renderer::add_line(const vec3f& start, const vec3f& end) 
{
  m_Vertices.emplace_back(start.x());
  m_Vertices.emplace_back(start.y());
  m_Vertices.emplace_back(start.z());
  m_Vertices.emplace_back(0.);
  m_Vertices.emplace_back(0.);
  m_Vertices.emplace_back(0.);

  m_Vertices.emplace_back(end.x());
  m_Vertices.emplace_back(end.y());
  m_Vertices.emplace_back(end.z());
  m_Vertices.emplace_back(0.);
  m_Vertices.emplace_back(0.);
  m_Vertices.emplace_back(0.);
}

void Line_renderer::add_line(const vec3f& start, const vec3f& end, const vec3f& color) 
{
  m_Vertices.emplace_back(start.x());
  m_Vertices.emplace_back(start.y());
  m_Vertices.emplace_back(start.z());
  m_Vertices.emplace_back(color.x());
  m_Vertices.emplace_back(color.y());
  m_Vertices.emplace_back(color.z());

  m_Vertices.emplace_back(end.x());
  m_Vertices.emplace_back(end.y());
  m_Vertices.emplace_back(end.z());
  m_Vertices.emplace_back(color.x());
  m_Vertices.emplace_back(color.y());
  m_Vertices.emplace_back(color.z());
}

void Line_renderer::draw()  
{
  assert(are_buffers_loaded());

  if (m_Vertices.empty()) return; 

  unsigned int nComponents = 6; // 6 components per vertex (3 for position + 3 for color)
  
  glBindVertexArray(m_VertexArray);

  glLineWidth(m_Width);
  glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(m_Vertices.size()/nComponents));

  glBindVertexArray(0);
  glLineWidth(1.f);
}

#endif // CGAL_GLFW_INTERNAL_LINE_RENDERER_H
