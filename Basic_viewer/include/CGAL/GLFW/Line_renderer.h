#ifndef LINE_RENDERER_H
#define LINE_RENDERER_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <vector>

class Line_renderer {
public:
  Line_renderer() : m_width(1.0f) {}

  ~Line_renderer() 
  {
    if (m_vao != 0) glDeleteVertexArrays(1, &m_vao);
    if (m_vbo[VBO_VERTEX] != 0) glDeleteBuffers(2, m_vbo);
  }

  inline 
  const std::vector<float>& get_vertices() const { return m_vertices; }

  inline 
  void set_width(const float width) { m_width = width; }

  inline 
  void initialize_buffers() 
  {
    glGenVertexArrays(1, &m_vao);
    glGenBuffers(2, m_vbo);

    glBindVertexArray(m_vao);

    glBindBuffer(GL_ARRAY_BUFFER, m_vbo[VBO_VERTEX]);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, m_vbo[VBO_COLOR]);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);  

    m_areBuffersInitialized = true;
  }

  inline 
  void load_buffers() 
  {
    assert(m_areBuffersInitialized);

    glBindVertexArray(m_vao);

    glBindBuffer(GL_ARRAY_BUFFER, m_vbo[VBO_VERTEX]);
    glBufferData(GL_ARRAY_BUFFER, m_vertices.size() * sizeof(float), m_vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, m_vbo[VBO_COLOR]);
    glBufferData(GL_ARRAY_BUFFER, m_colors.size() * sizeof(vec3f), m_colors.data(), GL_STATIC_DRAW);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    m_areBuffersLoaded = true;
  }

  inline  
  void add_line(const vec3f& start, const vec3f& end, const vec3f& color) 
  {
    m_colors.push_back(color);
    m_colors.push_back(color);

    m_vertices.push_back(start.x());
    m_vertices.push_back(start.y());
    m_vertices.push_back(start.z());

    m_vertices.push_back(end.x());
    m_vertices.push_back(end.y());
    m_vertices.push_back(end.z());
  }

  inline void draw()  
  {
    assert(m_areBuffersLoaded);

    if (m_vertices.empty()) return; 

    glBindVertexArray(m_vao);

    glLineWidth(m_width);
    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(m_vertices.size() / 3));

    glBindVertexArray(0);
    glLineWidth(1.f);
  }

private: 
  float m_width;

  std::vector<float> m_vertices;
  std::vector<vec3f> m_colors;

  enum VBOEnum
  {
    VBO_VERTEX=0, 
    VBO_COLOR, 
    NB_VBO
  };

  unsigned int m_vao, m_vbo[NB_VBO];

  bool m_areBuffersLoaded = false;
  bool m_areBuffersInitialized = false;
};

#endif // LINE_RENDERER_H



