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
class QOpenGLFramebufferObject;
namespace CGAL {
namespace Three {
//!
//! \brief The Primitive_container struct provides a base for the OpenGL data wrappers.
//!
struct DEMO_FRAMEWORK_EXPORT Primitive_container
{

  //!
  //! \brief Primitive_container constructor.
  //! \param program the `QOpenGLShaderProgram` used by the VAOs.
  //! \param indexed must be `true` if the data is indexed, `false` otherwise.
  //!
    Primitive_container(int program, bool indexed)
      : program_id(program), indexed(indexed),
        flat_size(0), idx_size(0)
    {}

    virtual ~Primitive_container()
    {
      Q_FOREACH(Vbo* vbo, VBOs)
        if(vbo)
          delete vbo;
      Q_FOREACH(CGAL::Three::Viewer_interface*viewer, VAOs.keys())
      {
        removeViewer(viewer);
      }
    }

    /*!
     * \brief bindUniformValues sets the uniform variables for the concerned shaders.
     *
     * Such variables are valid at every step of the pipeline. For example,
     * the ModelViewProjection matrix, the uniform color or the is_selected state are uniform values.
     * This function is called in the `draw()`function.
     * \attention `Vbo`s data should be allocated for this function to be effective.
     * \attention This should only be called once the `Vao`s and the `Vbo`s are created, in a
     * valid OpenGL context.
     * \param viewer the active `Viewer_interface`
     */
    void bindUniformValues(CGAL::Three::Viewer_interface* viewer) const
    {
      viewer->bindUniformValues(program_id);
      viewer->getShaderProgram(program_id)->bind();
      if(is_selected)
          viewer->getShaderProgram(program_id)->setUniformValue("is_selected", true);
      else
          viewer->getShaderProgram(program_id)->setUniformValue("is_selected", false);

      QColor c = color;
      if(program_id == Viewer_interface::PROGRAM_WITH_TEXTURE)
      {
         if(is_selected) c = c.lighter(120);
         viewer->getShaderProgram(program_id)->setAttributeValue
           ("color_facets",
            c.redF(),
            c.greenF(),
            c.blueF());
      }
      else if(program_id == Viewer_interface::PROGRAM_WITH_TEXTURED_EDGES)
      {
          if(is_selected) c = c.lighter(50);
          viewer->getShaderProgram(program_id)->setUniformValue
            ("color_lines",
             QVector3D(c.redF(), c.greenF(), c.blueF()));
      }
      viewer->getShaderProgram(program_id)->release();
    }

    //!
    //! \brief draw is the function that actually renders the data.
    //! \param viewer the active `Viewer_interface`.
    //! \param is_color_uniform should be `true` if the item is unicolor.
    //! \param fbo holds the texture that is used for transparency.
    //!
    virtual void draw(CGAL::Three::Viewer_interface* viewer,
                      bool is_color_uniform,
                      QOpenGLFramebufferObject* fbo = NULL) const = 0;

    //!
    //! \brief initializeBuffers sends the data to the GPU memory.
    //!
    //! It actually fills up the buffers with the data provided by `Vbo::allocate()`;
    //! \param viewer the active `Viewer_interface`.
    //!
    void initializeBuffers(CGAL::Three::Viewer_interface* viewer) const
    {
      if(!VAOs[viewer])
        return;
      viewer->makeCurrent();
      VAOs[viewer]->bind();
      Q_FOREACH(CGAL::Three::Vbo* vbo, VAOs[viewer]->vbos)
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
            VAOs[viewer]->program->enableAttributeArray(vbo->attribute);
            VAOs[viewer]->program->setAttributeBuffer(vbo->attribute, vbo->data_type, vbo->offset, vbo->tupleSize, vbo->stride);
          }
        }
        else if(vbo->vbo_type == QOpenGLBuffer::VertexBuffer)
        {
          VAOs[viewer]->program->disableAttributeArray(vbo->attribute);
        }
        vbo->release();
      }
      VAOs[viewer]->release();

      is_init[viewer] = true;
    }

    //!
    //! \brief initGL initializes the OpenGL containers.
    //! \attention It must be called within a valid OpenGL context. The `draw()` function of an item is always a safe place to call this.
    //!
    //! \param viewer the active `Viewer_interface`.
    virtual void initGL(CGAL::Three::Viewer_interface* viewer) const = 0;

    //!
    //! \brief removeViewer deletes and removes the Vao assigned to `viewer` from `Vaos`.
    //! \param viewer the `Viewer_interface` to remove.
    //!
    void removeViewer(CGAL::Three::Viewer_interface* viewer) const
    {
      delete VAOs[viewer];
      VAOs.remove(viewer);
    }

    //!
    //! \brief reset_vbos de-allocates the `Vbo`s. It must be called when the `Vbo`s data is updated.
    //!
    void reset_vbos()
    {
      Q_FOREACH(CGAL::Three::Vbo* vbo, VBOs)
      {
        if(!vbo)
          continue;
        vbo->allocated = false;
      }
    }

    //!
    //! \brief VAOs holds the `Vao`s for each `Viewer_interface`. As a `Vao` is context dependent, there must be one Vao for each `Viewer_interface`.
    //!
    mutable QMap<CGAL::Three::Viewer_interface*, Vao*> VAOs;
    //!
    //! \brief VBOs holds the `Vbo`s containing the data for this container.
    //!
    mutable std::vector<Vbo*> VBOs;
    //!
    //! \brief program_id is the `OpenGL_program_IDs` used with this container.
    //!
    int program_id;
    //!
    //! \brief indexed specifies if the data is indexed or not. This matters for the internal drawing functions.
    //!
    bool indexed;
    mutable QMap<CGAL::Three::Viewer_interface*, bool> is_init;
    mutable QMap<CGAL::Three::Viewer_interface*, bool> is_gl_init;

    //!
    //! \brief is_selected must be filled with the selected state of the item that holds this container everytime this state changes.
    //! If program_id doesn't use this property, it can be ignored.
    //!
    mutable bool is_selected;
    //!
    //! \brief flat_size contains the number of units contained in the vertex buffer.
    //! You can ignore it if `indexed` is true.
    //!
    mutable std::size_t flat_size;
    //! \brief flat_size contains the number of units contained in the barycenter buffer.
    //! You can ignore it if `program_id` is not an instanced program (like PROGRAM_SPHERES, PROGRAM_CUTPLANE_SPHERES,
    //! PROGRAM_INSTANCED_WIRE or PROGRAM_INSTANCED).
    mutable std::size_t center_size;
    //!
    //! \brief idx_size contains the number of indices in an `Index_buffer`. You can ignore it if `indexed` is `false`.
    //!
    mutable std::size_t idx_size;
    //!
    //! \brief color contains the color of the data. Ignore it if `is_color_uniform` is `false` in `draw()`.
    //!
    mutable QColor color;
}; //end of class Triangle_container

}
}

#endif // PRIMITIVE_CONTAINER_H
