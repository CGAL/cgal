// Copyright (c) 2017  GeometryFactory Sarl (France)
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
//
//
// Author(s)     : Maxime Gimeno

#ifndef PRIMITIVE_CONTAINER_H
#define PRIMITIVE_CONTAINER_H

#include <CGAL/license/Three.h>
#include <CGAL/Three/Buffer_objects.h>

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>

using namespace CGAL::Three;

#ifdef demo_framework_EXPORTS
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_EXPORT
#else
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_IMPORT
#endif
namespace CGAL {
namespace Three {

struct DEMO_FRAMEWORK_EXPORT Primitive_container
{

    Primitive_container(int program, bool indexed)
      : program_id(program), indexed(indexed),
        is_init(false), is_gl_init(false),
        flat_size(0), idx_size(0)
    {}

    virtual ~Primitive_container()
    {
      Q_FOREACH(Vbo* vbo, VBOs)
        if(vbo)
          delete vbo;
      delete VAO;
    }
    virtual void initGL(CGAL::Three::Scene_item* item, CGAL::Three::Viewer_interface* viewer) const = 0;

    virtual void draw(const Scene_item &item, CGAL::Three::Viewer_interface* viewer,
                      bool is_color_uniform) const = 0;

    void initializeBuffers() const
    {
      if(!VAO)
        return;
      VAO->bind();
      Q_FOREACH(CGAL::Three::Vbo* vbo, VAO->vbos)
      {
        vbo->bind();
        if(vbo->dataSize !=0)
        {
          if(!vbo->allocated)
          {
            if(vbo->vbo_type == QOpenGLBuffer::IndexBuffer)
              vbo->vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
            vbo->vbo.allocate(vbo->data, vbo->dataSize);
            vbo->allocated = true;
          }
          if(vbo->vbo_type == QOpenGLBuffer::VertexBuffer)
          {
            VAO->program->enableAttributeArray(vbo->attribute);
            VAO->program->setAttributeBuffer(vbo->attribute, vbo->data_type, vbo->offset, vbo->tupleSize, vbo->stride);
          }
        }
        else if(vbo->vbo_type == QOpenGLBuffer::VertexBuffer)
        {
          VAO->program->disableAttributeArray(vbo->attribute);
        }
        vbo->release();
      }
      VAO->release();

      is_init = true;
    }

    void reset_vbos()
    {
      Q_FOREACH(CGAL::Three::Vbo* vbo, VBOs)
      {
        if(!vbo)
          continue;
        vbo->allocated = false;
      }
    }

    mutable Vao* VAO;
    mutable std::vector<Vbo*> VBOs;
    int program_id;
    bool indexed;
    mutable bool is_init;
    mutable bool is_gl_init;
    mutable std::size_t flat_size;
    mutable std::size_t idx_size;
    mutable QColor color;
}; //end of class Triangle_container

}
}

#endif // PRIMITIVE_CONTAINER_H
