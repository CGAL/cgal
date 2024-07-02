#ifndef LINE_RENDERER_GLFW_H
#define LINE_RENDERER_GLFW_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <vector>

#include "utils.h"

class Line_renderer {
public:
  ~Line_renderer();

  void initialize_buffers();
  void load_buffers();
  void add_line(const vec3f& start, const vec3f& end, const vec3f& color);
  void draw();

  inline void set_width(const float width) { m_width = width; }

private: 
  std::vector<vec3f> m_data {};

  float m_width { 1.0f };

  unsigned int m_vao {};
  unsigned int m_vbo {};

  bool m_areBuffersLoaded { false };
  bool m_areBuffersInitialized { false };
};

/********************METHOD DEFINITIONS********************/

Line_renderer::~Line_renderer() 
{
  if (m_vao != 0) glDeleteVertexArrays(1, &m_vao);
  if (m_vbo != 0) glDeleteBuffers(1, &m_vbo);
}

inline 
void Line_renderer::initialize_buffers() 
{
  glGenVertexArrays(1, &m_vao);
  glGenBuffers(1, &m_vbo);

  glBindVertexArray(m_vao);

  glBindBuffer(GL_ARRAY_BUFFER, m_vbo);

  // vertex attribute  
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), nullptr);

  // color attribute 
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (const void*)12);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);  

  m_areBuffersInitialized = true;
}

inline 
void Line_renderer::load_buffers() 
{
  assert(m_areBuffersInitialized);

  glBindVertexArray(m_vao);

  glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
  glBufferData(GL_ARRAY_BUFFER, m_data.size() * sizeof(vec3f), m_data.data(), GL_STATIC_DRAW);

  glBindVertexArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  m_areBuffersLoaded = true;
}

inline  
void Line_renderer::add_line(const vec3f& start, const vec3f& end, const vec3f& color) 
{
  m_data.push_back(start);
  m_data.push_back(color);

  m_data.push_back(end);
  m_data.push_back(color);
}

inline 
void Line_renderer::draw()  
{
  assert(m_areBuffersLoaded);

  if (m_data.empty()) return; 

  glBindVertexArray(m_vao);

  glLineWidth(m_width);
  glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(m_data.size() / 2));

  glBindVertexArray(0);
  glLineWidth(1.f);
}

#endif // LINE_RENDERER_GLFW_H



