// Copyright (c) 2018  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

  //!
  //! \brief Creates a `Vao` from another one.
  //! \param program the `QOpenGLShaderProgram` corresponding to the one of `vao` but from
  //! the right viewer.
  //! \param vao the Vao to copy.
  //!
  //! All `vao`'s vbos will be shared. Use it for shared viewers.
  //! `initializeBuffers()` will have to be called again.
  //!
  //! \attention This must be called within a valid OpenGLContext.
  //! Most of the time, initGL() functions are safe places to do so.
  //!
  Vao( Vao* vao, QOpenGLShaderProgram* program)
    :vao(new QOpenGLVertexArrayObject()),
      program(program)
  {
    this->vao->create();
    vbos = vao->vbos;
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
    NORMALS,
    NOT_INSTANCED
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
