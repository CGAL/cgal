// Copyright (c) 2023, 2024  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include <iostream>

#include "Tools.h"
#include "Shader_program.h"

std::string Shader_program::s_shader_path;

//! \brief
void Shader_program::set_shader_path(const char* path)
{ s_shader_path = std::string(path); }

//! \brief
bool Shader_program::init() {
  initializeOpenGLFunctions();

  m_program = glCreateProgram();
  if (! m_program) {
    std::cout << "error creating shader program!\n";
    return false;
  }

  return true;
}

//! \brief
bool Shader_program::init(const char* vs, const char* gs, const char* fs) {
  if (! init()) return false;

  add_shader_from_file(vs, GL_VERTEX_SHADER);
  add_shader_from_file(gs, GL_GEOMETRY_SHADER);
  add_shader_from_file(fs, GL_FRAGMENT_SHADER);

  link();
  validate();

  return true;
}

//! \brief
bool Shader_program::init(const std::string& vs,
                          const std::string& gs,
                          const std::string& fs)
{ return init(vs.c_str(), gs.c_str(), fs.c_str()); }

//! \brief
bool Shader_program::init_with_vs_fs(const char* shader_file_prefix) {
  std::string prefix(shader_file_prefix);
  std::string vs = s_shader_path + prefix + "_vs.glsl";
  std::string fs = s_shader_path + prefix + "_fs.glsl";
  return init(vs, "", fs);
}

//! \brief
void Shader_program::add_shader(const char* shader_code, GLenum shader_type) {
  GLuint the_shader = glCreateShader(shader_type);

  const GLchar* the_code[] = { shader_code };
  GLint code_length[] = { static_cast<GLint>(strlen(shader_code)) };

  glShaderSource(the_shader, 1, the_code, code_length);
  glCompileShader(the_shader);


  GLint result = 0;
  GLchar elog[1024] = { 0 };
  glGetShaderiv(the_shader, GL_COMPILE_STATUS, &result);
  if (! result) {
    std::string shader_type_name;
    switch (shader_type) {
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

//! \brief
void Shader_program::add_shader_from_file(const char* shader_file,
                                          GLenum shader_type) {
  if (strlen(shader_file) == 0) return;

  auto src = read_file(shader_file);
  add_shader(src.c_str(), shader_type);
}

//! \brief
bool Shader_program::link() {
  GLint result = 0;
  GLchar elog[1024] = { 0 };

  glLinkProgram(m_program);
  glGetProgramiv(m_program, GL_LINK_STATUS, &result);
  if (! result) {
    glGetProgramInfoLog(m_program, sizeof(elog), NULL, elog);
    std::cout << "! error linking program:\n" << elog << std::endl;
    return false;
  }
  return true;
}

//! \brief
bool Shader_program::validate() {
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

//! \brief
GLint Shader_program::get_uniform_location(const GLchar* name) {
  const auto uniform_loc = glGetUniformLocation(m_program, name);
  return uniform_loc;
}

//! \brief
void Shader_program::use() { glUseProgram(m_program); }

//! \brief
void Shader_program::unuse() { glUseProgram(0); }

//! \brief
void Shader_program::set_uniform(GLint uniform_loc, const QMatrix4x4& m)
{ glUniformMatrix4fv(uniform_loc, 1, GL_FALSE, m.data()); }

//! \brief
void Shader_program::set_uniform(const GLchar* name, const QMatrix4x4& m) {
  const auto uniform_loc = get_uniform_location(name);
  set_uniform(uniform_loc, m);
}

//! \brief
void Shader_program::set_uniform(GLint uniform_loc, const QMatrix3x3& m)
{ glUniformMatrix3fv(uniform_loc, 1, GL_FALSE, m.data()); }

//! \brief
void Shader_program::set_uniform(const GLchar* name, const QMatrix3x3& m) {
  const auto uniform_loc = get_uniform_location(name);
  set_uniform(uniform_loc, m);
}

//! \brief
void Shader_program::set_uniform(GLint uniform_loc, const QVector4D& v)
{ glUniform4f(uniform_loc, v.x(), v.y(), v.z(), v.w()); }

//! \brief
void Shader_program::set_uniform(const GLchar* name, const QVector4D& v) {
  const auto uniform_loc = get_uniform_location(name);
  set_uniform(uniform_loc, v);
}
