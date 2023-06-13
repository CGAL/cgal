
#include "Shader_program.h"

#include <iostream>

#include "Tools.h"


bool Shader_program::init()
{
  initializeOpenGLFunctions();

  m_program = glCreateProgram();
  if (!m_program)
  {
    std::cout << "error creating shader program!\n";
    return false;
  }

  return true;
}
bool Shader_program::init(const char* vs, const char* gs, const char* fs)
{
  return false;
}
void Shader_program::add_shader(const char* shader_code, GLenum shader_type)
{
  GLuint the_shader = glCreateShader(shader_type);

  const GLchar* the_code[] = { shader_code };
  GLint code_length[] = { strlen(shader_code) };

  glShaderSource(the_shader, 1, the_code, code_length);
  glCompileShader(the_shader);


  GLint result = 0;
  GLchar elog[1024] = { 0 };
  glGetShaderiv(the_shader, GL_COMPILE_STATUS, &result);
  if (!result)
  {
    std::string shader_type_name;
    switch (shader_type)
    {
    case GL_VERTEX_SHADER:   shader_type_name = "VERTEX"; break;
    case GL_GEOMETRY_SHADER: shader_type_name = "GEOMETRY"; break;
    case GL_FRAGMENT_SHADER: shader_type_name = "FRAGMENT"; break;
    }
    glGetShaderInfoLog(the_shader, sizeof(elog), NULL, elog);
    std::cout << "! error compiling the " << shader_type_name <<
      " shader:\n" << elog << std::endl;
    return;
  }

  glAttachShader(m_program, the_shader);
}
void Shader_program::add_shader_from_file(const char* shader_file, 
                                          GLenum shader_type)
{
  auto src = read_file(shader_file);
  add_shader(src.c_str(), shader_type);
}

bool Shader_program::link()
{
  GLint result = 0;
  GLchar elog[1024] = { 0 };

  glLinkProgram(m_program);
  glGetProgramiv(m_program, GL_LINK_STATUS, &result);
  if (!result)
  {
    glGetProgramInfoLog(m_program, sizeof(elog), NULL, elog);
    std::cout << "! error linking program:\n" << elog << std::endl;
    return false;
  }
  return true;
}
bool Shader_program::validate()
{
  GLint result = 0;
  GLchar elog[1024] = { 0 };

  glValidateProgram(m_program);
  glGetProgramiv(m_program, GL_VALIDATE_STATUS, &result);
  if (!result)
  {
    glGetProgramInfoLog(m_program, sizeof(elog), NULL, elog);
    std::cout << "! error validating program:\n" << elog << std::endl;
    return false;
  }
  return true;
}

GLint Shader_program::get_uniform_location(const GLchar* name)
{
  const auto uniform_loc = glGetUniformLocation(m_program, name);
  return uniform_loc;
}

void Shader_program::use()
{
  glUseProgram(m_program);
}
void Shader_program::unuse()
{
  glUseProgram(0);
}

void Shader_program::set_uniform(GLint uniform_loc, const QMatrix4x4& m)
{
  glUniformMatrix4fv(uniform_loc, 1, GL_FALSE, m.data());
}
void Shader_program::set_uniform(const GLchar* name, const QMatrix4x4& m)
{
  const auto uniform_loc = get_uniform_location(name);
  set_uniform(uniform_loc, m);
}

void Shader_program::set_uniform(GLint uniform_loc, const QVector4D& v)
{
  glUniform4f(uniform_loc, v.x(), v.y(), v.z(), v.w());
}
void Shader_program::set_uniform(const GLchar* name, const QVector4D& v)
{
  const auto uniform_loc = get_uniform_location(name);
  set_uniform(uniform_loc, v);
}