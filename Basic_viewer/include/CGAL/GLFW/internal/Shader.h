#ifndef CGAL_SHADER_GLFW_H
#define CGAL_SHADER_GLFW_H

#include <unordered_map>
#include <iostream>
#include <memory>

#include <glad/glad.h>

class Shader
{
public:
  Shader()=default;
  Shader(int program);
  ~Shader();

  inline void use() const { glUseProgram(m_Program); }
  inline void destroy() { if (m_Program != 0) glDeleteProgram(m_Program); }

  int get_uniform_location(const GLchar* name);

  inline void set_mat4f(const GLchar* name, const GLfloat* data, GLboolean transpose=false) { glUniformMatrix4fv(get_uniform_location(name), 1, transpose, data); }
  inline void set_vec4f(const GLchar* name, const GLfloat* data) { glUniform4fv(get_uniform_location(name), 1, data); }
  inline void set_vec3f(const GLchar* name, const GLfloat* data) { glUniform3fv(get_uniform_location(name), 1, data); }
  inline void set_vec2f(const GLchar* name, const GLfloat* data) { glUniform2fv(get_uniform_location(name), 1, data); }
  inline void set_float(const GLchar* name, const float data) { glUniform1f(get_uniform_location(name), data); }
  inline void set_int(const GLchar* name, const int data) { glUniform1i(get_uniform_location(name), data); }

  static std::shared_ptr<Shader> create(const GLchar* sourceVertex, const GLchar* sourceFragment, const GLchar* sourceGeometry=nullptr);

  static GLchar* enum_type_to_str(GLenum type);

private:
  static void check_compile_errors(unsigned int  shader, const GLchar* type);
  static unsigned int  compile_shader(GLenum type, const GLchar* source);

private:
  std::unordered_map<const GLchar*, int> m_Uniforms {};
  int m_Program { 0 };
};

/********************METHOD IMPLEMENTATIONS********************/

Shader::Shader(int program) : 
  m_Program(program) 
{
}

Shader::~Shader()
{
  destroy();
}

int Shader::get_uniform_location(const GLchar* name)
{
  int loc = m_Uniforms[name];
  if (loc != 0)
  {
    return loc;
  }

  loc = glGetUniformLocation(m_Program, name);
  m_Uniforms[name] = loc;
  return loc;
}

GLchar* Shader::enum_type_to_str(GLenum type)
{
  switch (type)
  {
  case GL_VERTEX_SHADER:   return "VERTEX";
  case GL_GEOMETRY_SHADER: return "GEOMETRY";
  case GL_FRAGMENT_SHADER: return "FRAGMENT";
  default: 
    std::cerr << "Unsupported shader enum type!" << std::endl; 
    return "UNSUPPORTED TYPE";
  }
}

unsigned int  Shader::compile_shader(GLenum type, const GLchar* source)
{
  unsigned int  shader = glCreateShader(type);
  glShaderSource(shader, 1, &source, nullptr);
  glCompileShader(shader);
  Shader::check_compile_errors(shader, Shader::enum_type_to_str(type));
  return shader; 
}

std::shared_ptr<Shader> Shader::create(const GLchar* sourceVertex, const GLchar* sourceFragment, const GLchar* sourceGeometry)
{
  unsigned int  vertexShader = Shader::compile_shader(GL_VERTEX_SHADER, sourceVertex);
  unsigned int  geometryShader;
  if (sourceGeometry)
  {
    geometryShader = Shader::compile_shader(GL_GEOMETRY_SHADER, sourceGeometry);
  }
  unsigned int  fragmentShader = Shader::compile_shader(GL_FRAGMENT_SHADER, sourceFragment);

  unsigned int  program = glCreateProgram();
  glAttachShader(program, vertexShader);
  if (sourceGeometry)
  {
    glAttachShader(program, geometryShader);
  }
  glAttachShader(program, fragmentShader);

  glLinkProgram(program);
  Shader::check_compile_errors(program, "PROGRAM");

  glDeleteShader(vertexShader);
  if (sourceGeometry)   
  {
    glDeleteShader(geometryShader);
  }
  glDeleteShader(fragmentShader);

  return std::make_shared<Shader>(program);
}

void Shader::check_compile_errors(unsigned int  shader, const GLchar* type)
{
  GLint success;
  GLchar infoLog[1024];

  if (type != "PROGRAM")
  {
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

    if (!success)
    {
      glGetShaderInfoLog(shader, 1024, NULL, infoLog);
      std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << "\n"
                << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
    }
  }
  else
  {
    glGetProgramiv(shader, GL_LINK_STATUS, &success);

    if (!success)
    {
      glGetProgramInfoLog(shader, 1024, NULL, infoLog);
      std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << "\n"
                << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
    }
  }
}

#endif // CGAL_SHADER_GLFW_H