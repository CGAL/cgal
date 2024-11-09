#ifndef CGAL_GLFW_INTERNAL_BUFFER_VAO_H
#define CGAL_GLFW_INTERNAL_BUFFER_VAO_H

#include <memory>

#include <glad/glad.h>

#include "BufferObject.h"

// Inspired by https://github.com/TheCherno/Hazel/blob/master/Hazel/src/Platform/OpenGL/OpenGLVertexArray.h
class VAO 
{
public:
  VAO();
  ~VAO();

  void set_VBO(const std::shared_ptr<VBO>& vbo);
  void set_EBO(const std::shared_ptr<EBO>& ebo);

  inline const std::shared_ptr<VBO>& get_VBO() const { return m_VBO; }
  inline const std::shared_ptr<EBO>& get_EBO() const { return m_EBO; }

  inline void bind()   const { glBindVertexArray(m_Id); }
  inline void unbind() const { glBindVertexArray(m_Id); }
  inline void destroy() {glDeleteVertexArrays(1, &m_Id); }

  inline static std::shared_ptr<VAO> create() { return std::make_shared<VAO>(); } 

private:
  unsigned int m_Id;
  unsigned int m_VBOAttribIndex;
  std::shared_ptr<VBO> m_VBO;
  std::shared_ptr<EBO> m_EBO; 
};

VAO::VAO() : m_VBOAttribIndex(0) 
{ 
  glGenVertexArrays(1, &m_Id); 
}

VAO::~VAO()
{
  destroy();
}

void VAO::set_VBO(const std::shared_ptr<VBO>& vbo) 
{
  glBindVertexArray(m_Id);
  vbo->bind();

  const auto& layout = vbo->get_layout(); 
  for (const auto& element : layout)
  {
    switch (element.Type)
    {
      case ShaderDataType::FLOAT:
      case ShaderDataType::FLOAT2:
      case ShaderDataType::FLOAT3:
      case ShaderDataType::FLOAT4:
        glEnableVertexAttribArray(m_VBOAttribIndex);
        glVertexAttribPointer(m_VBOAttribIndex, 
                              element.get_component_count(), 
                              get_shader_data_type_to_opengl_base_type(element.Type),
                              element.Normalized ? GL_TRUE : GL_FALSE, 
                              layout.get_stride(), 
                              (const void*) element.Offset);
        m_VBOAttribIndex++;
        break;
      case ShaderDataType::INT:
      case ShaderDataType::INT2:
      case ShaderDataType::INT3:
      case ShaderDataType::INT4:
      case ShaderDataType::UINT:
      case ShaderDataType::UINT2:
      case ShaderDataType::UINT3:
      case ShaderDataType::UINT4:
      case ShaderDataType::BOOL:
        glEnableVertexAttribArray(m_VBOAttribIndex);
        glVertexAttribIPointer(m_VBOAttribIndex, 
                                element.get_component_count(), 
                                get_shader_data_type_to_opengl_base_type(element.Type),
                                layout.get_stride(), 
                                (const void*) element.Offset);
        m_VBOAttribIndex++;
        break;
      case ShaderDataType::MAT3:
      case ShaderDataType::MAT4:
        unsigned int count = element.get_component_count();
        for (unsigned int i = 0; i < count; ++i)
        {
          glEnableVertexAttribArray(m_VBOAttribIndex);
          glVertexAttribPointer(m_VBOAttribIndex, 
                                count, 
                                get_shader_data_type_to_opengl_base_type(element.Type), 
                                element.Normalized ? GL_TRUE : GL_FALSE, 
                                layout.get_stride(), 
                                (const void*) (element.Offset + sizeof(float) * count * i));
          glVertexAttribDivisor(m_VBOAttribIndex, 1);
          m_VBOAttribIndex++;
        }
        break;
    }
  }
  m_VBO = vbo;
}

void VAO::set_EBO(const std::shared_ptr<EBO>& ebo) 
{ 
  glBindVertexArray(m_Id);
  ebo->bind(); 

  m_EBO = ebo; 
}

#endif // CGAL_GLFW_INTERNAL_BUFFER_VAO_H
