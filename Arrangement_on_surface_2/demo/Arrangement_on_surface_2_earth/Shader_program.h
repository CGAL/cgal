// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef SHADER_PROGRAM_H
#define SHADER_PROGRAM_H

#include <string>
#include <qmatrix4x4.h>

#include "Common_defs.h"

class Shader_program : protected OpenGLFunctionsBase {
public:
  static void set_shader_path(const char* path);

  bool init();
  bool init(const char* vs, const char* gs, const char* fs);
  bool init(const std::string& vs, const std::string& gs,
            const std::string& fs);

  /*! Initialize with just the vertex and fragment shaders
   */
  bool init_with_vs_fs(const char* shader_file_prefix);

  void add_shader(const char* shader_code, GLenum shader_type);
  void add_shader_from_file(const char* shader_file, GLenum shader_type);

  bool link();
  bool validate();

  GLint get_uniform_location(const GLchar* name);

  void use();
  void unuse();

  void set_uniform(GLint uniform_loc, const QMatrix4x4& m);
  void set_uniform(const GLchar* name, const QMatrix4x4& m);

  void set_uniform(GLint uniform_loc, const QMatrix3x3& m);
  void set_uniform(const GLchar* name, const QMatrix3x3& m);

  void set_uniform(GLint uniform_loc, const QVector4D& v);
  void set_uniform(const GLchar* name, const QVector4D& v);

private:
  GLuint m_program;
  static std::string s_shader_path;
};


#endif
