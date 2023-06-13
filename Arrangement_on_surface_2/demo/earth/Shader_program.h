
#ifndef SHADER_PROGRAM_H
#define SHADER_PROGRAM_H

#include <qmatrix4x4.h>
#include "Common_defs.h"


class Shader_program : protected OpenGLFunctionsBase
{
public:

  bool init();
  bool init(const char* vs, const char* gs, const char* fs);

  void add_shader(const char* shader_code, GLenum shader_type);
  void add_shader_from_file(const char* shader_file, GLenum shader_type);

  bool link();
  bool validate();

  GLint get_uniform_location(const GLchar* name);

  void use();
  void unuse();

  void set_uniform(GLint uniform_loc, const QMatrix4x4& m);
  void set_uniform(const GLchar* name, const QMatrix4x4& m);

  void set_uniform(GLint uniform_loc, const QVector4D& v);
  void set_uniform(const GLchar* name, const QVector4D& v);

//private:
  GLuint  m_program;
};


#endif
