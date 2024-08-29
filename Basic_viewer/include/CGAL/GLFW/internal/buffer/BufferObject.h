#ifndef CGAL_BUFFER_H
#define CGAL_BUFFER_H

#include <iostream>
#include <vector>
#include <memory>

#include <glad/glad.h>

// Inspired by https://github.com/TheCherno/Hazel/blob/master/Hazel/src/Platform/OpenGL/OpenGLBuffer.h
enum class ShaderDataType 
{
  FLOAT, FLOAT2, FLOAT3, FLOAT4, INT, INT2, INT3, INT4, UINT, UINT2, UINT3, UINT4, BOOL, MAT3, MAT4 
};

GLenum get_shader_data_type_to_opengl_base_type(ShaderDataType type)
{
  switch(type)
  {
    case ShaderDataType::BOOL:   return GL_BOOL;
    case ShaderDataType::FLOAT:  return GL_FLOAT;
    case ShaderDataType::FLOAT2: return GL_FLOAT;
    case ShaderDataType::FLOAT3: return GL_FLOAT;
    case ShaderDataType::FLOAT4: return GL_FLOAT;
    case ShaderDataType::MAT3:   return GL_FLOAT;
    case ShaderDataType::MAT4:   return GL_FLOAT;
    case ShaderDataType::INT:    return GL_INT;
    case ShaderDataType::INT2:   return GL_INT;
    case ShaderDataType::INT3:   return GL_INT;
    case ShaderDataType::INT4:   return GL_INT;
    case ShaderDataType::UINT:   return GL_UNSIGNED_INT;
    case ShaderDataType::UINT2:  return GL_UNSIGNED_INT;
    case ShaderDataType::UINT3:  return GL_UNSIGNED_INT;
    case ShaderDataType::UINT4:  return GL_UNSIGNED_INT;
    default: throw(errno);
  }
}

size_t get_shader_data_type_size(ShaderDataType type)
{
  switch(type)
  {
    case ShaderDataType::BOOL:   return 1;
    case ShaderDataType::FLOAT:  return 4;
    case ShaderDataType::FLOAT2: return 4 * 2;
    case ShaderDataType::FLOAT3: return 4 * 3;
    case ShaderDataType::FLOAT4: return 4 * 4;
    case ShaderDataType::INT:    return 4;
    case ShaderDataType::INT2:   return 4 * 2;
    case ShaderDataType::INT3:   return 4 * 3;
    case ShaderDataType::INT4:   return 4 * 4;
    case ShaderDataType::UINT:   return 4;
    case ShaderDataType::UINT2:  return 4 * 2;
    case ShaderDataType::UINT3:  return 4 * 3;
    case ShaderDataType::UINT4:  return 4 * 4;
    case ShaderDataType::MAT3:   return 4 * 3 * 3;
    case ShaderDataType::MAT4:   return 4 * 4 * 4;
    default: throw(errno);
  }
}

struct BufferElement
{
  ShaderDataType Type;
  size_t Size; 
  size_t Offset;
  bool Normalized; 

  BufferElement(ShaderDataType type, bool normalized=false) : Type(type), Size(get_shader_data_type_size(type)), Offset(0), Normalized(normalized) {}  

  unsigned int get_component_count() const
  {
    switch(Type)
    {
      case ShaderDataType::BOOL:   return 1;
      case ShaderDataType::FLOAT:  return 1;
      case ShaderDataType::FLOAT2: return 2;
      case ShaderDataType::FLOAT3: return 3;
      case ShaderDataType::FLOAT4: return 4;
      case ShaderDataType::INT:    return 1;
      case ShaderDataType::INT2:   return 2;
      case ShaderDataType::INT3:   return 3;
      case ShaderDataType::INT4:   return 4;
      case ShaderDataType::UINT:   return 1;
      case ShaderDataType::UINT2:  return 2;
      case ShaderDataType::UINT3:  return 3;
      case ShaderDataType::UINT4:  return 4;
      case ShaderDataType::MAT3:   return 3; // 3 * float3;
      case ShaderDataType::MAT4:   return 4; // 4 * float4;
      default: throw(errno);
    }
  }
};

class BufferLayout
{
public:
  BufferLayout() : m_Elements(), m_Stride(0) {} 
  BufferLayout(std::initializer_list<BufferElement> elements) : m_Elements(elements), m_Stride(0) 
  {
    compute_stride_and_offset();
  }

  inline size_t get_stride() const { return m_Stride; }

  inline std::vector<BufferElement>::iterator begin() { return m_Elements.begin(); }
  inline std::vector<BufferElement>::iterator end() { return m_Elements.end(); }
  inline std::vector<BufferElement>::const_iterator begin() const { return m_Elements.begin(); }
  inline std::vector<BufferElement>::const_iterator end() const { return m_Elements.end(); }

private: 
  void compute_stride_and_offset();

private: 
  std::vector<BufferElement> m_Elements; 
  size_t m_Stride; 
}; 

void BufferLayout::compute_stride_and_offset()
{
  size_t offset = 0;
  for (auto& element : m_Elements)
  {
    element.Offset = offset; 
    offset += element.Size; 
    m_Stride += element.Size; 
  }
}

/*********************************Vertex Buffer Object*********************************/

class VBO
{
public:
  VBO();
  VBO(float* vertices, unsigned int size);
  ~VBO();

  inline void load(float* vertices, unsigned int size)
  {
    glBindBuffer(GL_ARRAY_BUFFER, m_Id);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * size, vertices, GL_STATIC_DRAW);   
  }

  inline const BufferLayout& get_layout() const { return m_Layout; }
  inline void set_layout(const BufferLayout& layout) { m_Layout = layout; }
  
  inline void bind()   const { glBindBuffer(GL_ARRAY_BUFFER, m_Id); }
  inline void unbind() const { glBindBuffer(GL_ARRAY_BUFFER, 0); } 
  inline void destroy() { glDeleteBuffers(1, &m_Id); }

  inline static std::shared_ptr<VBO> create() { return std::make_shared<VBO>(); }
  inline static std::shared_ptr<VBO> create(float* vertices, unsigned int size) { return std::make_shared<VBO>(vertices, size); }

private:
  unsigned int m_Id;
  BufferLayout m_Layout;
};

VBO::VBO()
{ 
  glGenBuffers(1, &m_Id); 
}

VBO::VBO(float* vertices, unsigned int size)
{ 
  glGenBuffers(1, &m_Id); 
  glBindBuffer(GL_ARRAY_BUFFER, m_Id);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float) * size, vertices, GL_STATIC_DRAW);  
}

VBO::~VBO()
{ 
  destroy();
}

/*********************************Index Buffer Object*********************************/

class EBO
{
public:
  EBO();
  EBO(unsigned int* indices, unsigned int count);
  ~EBO();

  inline void load(unsigned int* indices, unsigned int count)
  {
    m_Count = count;
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_Id);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * m_Count, indices, GL_STATIC_DRAW);   
  }

  inline unsigned int get_count() const { return m_Count; }

  inline void bind()   const { glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_Id); }
  inline void unbind() const { glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); } 
  inline void destroy() { glDeleteBuffers(1, &m_Id); }

  inline static std::shared_ptr<EBO> create() { return std::make_shared<EBO>(); }
  inline static std::shared_ptr<EBO> create(unsigned int* indices, unsigned int count) { return std::make_shared<EBO>(indices, count); }

private:
  unsigned int m_Id;
  unsigned int m_Count; 
};

EBO::EBO()  
{ 
  glGenBuffers(1, &m_Id); 
}

EBO::EBO(unsigned int* indices, unsigned int count)  
{ 
  glGenBuffers(1, &m_Id); 
  m_Count = count;
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_Id);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * m_Count, indices, GL_STATIC_DRAW);   
}

EBO::~EBO()
{
  destroy();
}

#endif // CGAL_BUFFER_H
