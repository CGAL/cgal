#ifndef CGAL_SHADER_GLFW_H
#define CGAL_SHADER_GLFW_H

#include <unordered_map>
#include <iostream>

#include <glad/glad.h>

class Shader
{
public:
  Shader()=default;
  Shader(int program);

  void use();
  void destroy();

  int get_uniform_location(const GLchar* name);

  void set_mat4f(const GLchar* name, const GLfloat* data, GLboolean transpose=false);
  void set_vec4f(const GLchar* name, const GLfloat* data);
  void set_vec3f(const GLchar* name, const GLfloat* data);
  void set_float(const GLchar* name, const float data);
  void set_int(const GLchar* name, const int data);

  static Shader create_shader(const GLchar* sourceVertex, const GLchar* sourceFragment, const GLchar* sourceGeometry=nullptr);

  static GLchar* enum_type_to_str(GLenum type);

private:
  static void check_compile_errors(GLuint shader, const GLchar* type);
  static GLuint compile_shader(GLenum type, const GLchar* source);

private:
  std::unordered_map<const GLchar*, int> m_uniforms {};
  int m_program { 0 };
};

/********************METHOD DEFINITIONS********************/

Shader::Shader(int program) : 
  m_program(program) 
{
}

inline 
void Shader::use()
{
  glUseProgram(m_program);
}

inline 
void Shader::destroy()
{
  if (m_program != 0) glDeleteProgram(m_program);
}

inline 
int Shader::get_uniform_location(const GLchar* name)
{
  int loc = m_uniforms[name];
  if (loc != 0)
  {
    return loc;
  }

  loc = glGetUniformLocation(m_program, name);
  m_uniforms[name] = loc;
  return loc;
}

inline 
void Shader::set_mat4f(const GLchar* name, const GLfloat* data, GLboolean transpose)
{
  glUniformMatrix4fv(get_uniform_location(name), 1, transpose, data);
}

inline 
void Shader::set_vec4f(const GLchar* name, const GLfloat* data)
{
  glUniform4fv(get_uniform_location(name), 1, data);
}

inline 
void Shader::set_vec3f(const GLchar* name, const GLfloat* data)
{
  glUniform3fv(get_uniform_location(name), 1, data);
}

inline 
void Shader::set_float(const GLchar* name, const float data)
{
  glUniform1f(get_uniform_location(name), data);
}

inline
void Shader::set_int(const GLchar* name, const int data)
{
  glUniform1i(get_uniform_location(name), data);
}

inline 
GLchar* Shader::enum_type_to_str(GLenum type)
{
  switch (type)
  {
  case GL_VERTEX_SHADER:   return "VERTEX";
  case GL_GEOMETRY_SHADER: return "GEOMETRY";
  case GL_FRAGMENT_SHADER: return "FRAGMENT";
  default:                 return "UNSUPPORTED TYPE";
  }
}

inline 
GLuint Shader::compile_shader(GLenum type, const GLchar* source)
{
  GLuint shader = glCreateShader(type);
  glShaderSource(shader, 1, &source, nullptr);
  glCompileShader(shader);
  Shader::check_compile_errors(shader, Shader::enum_type_to_str(type));
  return shader; 
}

inline 
Shader Shader::create_shader(const GLchar* sourceVertex, const GLchar* sourceFragment, const GLchar* sourceGeometry)
{
  GLuint vertexShader = Shader::compile_shader(GL_VERTEX_SHADER, sourceVertex);
  GLuint geometryShader;
  if (sourceGeometry)
  {
    geometryShader = Shader::compile_shader(GL_GEOMETRY_SHADER, sourceGeometry);
  }
  GLuint fragmentShader = Shader::compile_shader(GL_FRAGMENT_SHADER, sourceFragment);

  GLuint program = glCreateProgram();
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

  return Shader(program);
}

inline 
void Shader::check_compile_errors(GLuint shader, const GLchar* type)
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