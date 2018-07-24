// Copyright (c) 2018  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Maxime Gimeno

#ifndef BUFFER_OBJECTS_H
#define BUFFER_OBJECTS_H

#include <CGAL/license/Three.h>


#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>

#include <vector>

namespace CGAL {
namespace Three {

struct Vbo;
//!
//! \brief The Vao struct is a wrapper for the `QOpenGLVertexArrayObject`.
//! They are context dependent, most of the time it means `Viewer` dependent.
//!
struct Vao{

  QOpenGLVertexArrayObject* vao;
  QOpenGLShaderProgram* program;
  std::vector<Vbo*> vbos;
  //!
  //! \brief Creates a `Vao`.
  //! \param program the `QOpenGLShaderProgram` associated with this `Vao`.
  //! \attention This must be called within a valid OpenGLContext.
  //! Most of the time, initGL() functions are safe places to do so.
  //!
  Vao( QOpenGLShaderProgram* program)
    :vao(new QOpenGLVertexArrayObject()),
      program(program)
  {
    vao->create();
  }

  ~Vao()
  {
    delete vao;
  }

  //!Adds a `Vbo` to the list of Vbos.
  void addVbo(Vbo* vbo)
  {
    vbos.push_back(vbo);
  }

  //!Makes the `Vao` and its program active until `release()` is called.
  void bind()
  {
    program->bind();
    vao->bind();
  }

  //!Makes the `Vao` and its program not active.
  void release()
  {
    vao->release();
    program->release();
  }
};

//!
//! \brief The Vbo struct is a wrapper for the `QOpenGLBufferObject` item.
//! It contains the data necessary for the rendering of any displayed entity.
//! A Vbo can be shared between Vaos of the same context.
struct Vbo
{
enum Flag{
  GEOMETRY=0,
  COLORS,
  NORMALS
};
  QOpenGLBuffer vbo;
  const char* attribute;
  Flag flag;
  void* data;
  int dataSize;
  QOpenGLBuffer::Type vbo_type;
  GLenum data_type;
  int offset;
  int tupleSize;
  int stride;
  bool allocated;
  //!
  //! \brief Creates a `Vbo`.
  //! \param attribute the name of the corresponding data in the shader.
  //! \param flag the flag that specifies which type of data this corresponds to.
  //! \param vbo_type is almost always `QOpenGLBuffer::VertexBuffer` but can be `QOpenGLBuffer::IndexBuffer`
  //! if it contains the indices for an indexed rendering.
  //! \param data_type the GL data type. Mostly GL_FLOAT.
  //! \param offset the offset in the buffer. See OpenGL documentation.
  //! \param tupleSize the size of the tuple. If it contains vector_3s, then tuple_size is 3.
  //! \param stride the stride for the buffer. See OpenGL documentation.

  //! \attention This must be called within a valid OpenGLContext.
  //! Most of the time, initGL() functions are safe places to do so.
  //!
  Vbo(
      const char* attribute,
      Flag flag,
      QOpenGLBuffer::Type vbo_type = QOpenGLBuffer::VertexBuffer,
      GLenum data_type=GL_FLOAT,
      int offset = 0,
      int tupleSize = 3,
      int stride = 0):
    attribute(attribute),
    flag(flag),
    data(0),
    dataSize(0),
    vbo_type(vbo_type),
    data_type(data_type),
    offset(offset),
    tupleSize(tupleSize),
    stride(stride),
    allocated(false)
  {
    vbo = QOpenGLBuffer(vbo_type);
    vbo.create();
  }
  ~Vbo()
  {
    vbo.destroy();
  }

  //! Makes the `Vbo` active until `release()` is called.
  bool bind()
  {
    return vbo.bind();
  }

  //!Makes the `Vbo` and its program not active.
  void release()
  {
    vbo.release();
  }

  //!
  //! \brief allocate gives CPU data to the GPU.
  //! \param data the content of the buffer.
  //! \param datasize the number of bytes in `data`.
  //!
  void allocate(void* data, int datasize)
  {
    this->data = data;
    this->dataSize = datasize;
  }

};

}
}
#endif // BUFFER_OBJECTS_H
