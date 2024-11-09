#ifndef CGAL_GLFW_INTERNAL_LINE_RENDERER_H
#define CGAL_GLFW_INTERNAL_LINE_RENDERER_H

#include <vector>

#include <glad/glad.h>

#include "utils.h"
#include "buffer/VAO.h"

class Line_renderer {
public:
  void initialize_buffers(const BufferLayout& layout);
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

  std::shared_ptr<VAO> m_VAO;
  std::shared_ptr<VBO> m_VBO;

  bool m_AreBuffersLoaded { false };
  bool m_AreBuffersInitialized { false };
};

/********************METHOD IMPLEMENTATIONS********************/

void Line_renderer::initialize_buffers(const BufferLayout& layout) 
{
  m_VAO = VAO::create();
  m_VBO = VBO::create();
  m_VBO->set_layout(layout);
  m_VAO->set_VBO(m_VBO);

  m_AreBuffersInitialized = true;
}

void Line_renderer::load_buffers() 
{
  assert(are_buffers_initialized());

  m_VAO->bind();
  m_VBO->load(m_Vertices.data(), m_Vertices.size());

  m_AreBuffersLoaded = true;
}

void Line_renderer::add_line(const vec3f& start, const vec3f& end) 
{
  m_Vertices.push_back(start.x());
  m_Vertices.push_back(start.y());
  m_Vertices.push_back(start.z());

  m_Vertices.push_back(end.x());
  m_Vertices.push_back(end.y());
  m_Vertices.push_back(end.z());
}

void Line_renderer::add_line(const vec3f& start, const vec3f& end, const vec3f& color) 
{
  m_Vertices.push_back(start.x());
  m_Vertices.push_back(start.y());
  m_Vertices.push_back(start.z());
  m_Vertices.push_back(color.x());
  m_Vertices.push_back(color.y());
  m_Vertices.push_back(color.z());

  m_Vertices.push_back(end.x());
  m_Vertices.push_back(end.y());
  m_Vertices.push_back(end.z());
  m_Vertices.push_back(color.x());
  m_Vertices.push_back(color.y());
  m_Vertices.push_back(color.z());
}

void Line_renderer::draw()  
{
  assert(are_buffers_loaded());

  if (m_Vertices.empty()) return; 

  unsigned int nComponents = m_VBO->get_layout().get_stride() / sizeof(float);
  
  m_VAO->bind();
  glLineWidth(m_Width);
  glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(m_Vertices.size()/nComponents)); // 6 = n components per vertex
  glLineWidth(1.f);
}

#endif // CGAL_GLFW_INTERNAL_LINE_RENDERER_H



