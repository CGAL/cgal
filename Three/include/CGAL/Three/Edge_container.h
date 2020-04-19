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

#ifndef EDGE_CONTAINER_H
#define EDGE_CONTAINER_H

#include <CGAL/license/Three.h>


#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Primitive_container.h>

using namespace CGAL::Three;
#ifdef demo_framework_EXPORTS
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_EXPORT
#else
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_IMPORT
#endif
struct Edge_d;
namespace CGAL {
namespace Three {

//!
//! \brief The Edge_container struct wraps the OpenGL data for drawing lines.
//!
struct DEMO_FRAMEWORK_EXPORT Edge_container :public Primitive_container
{

  //!
  //! \brief The vbosName enum
  //!
  //! Holds the `Vbo` Ids of this container.
  //!
  enum vbosName {
    Vertices = 0, //!< Designates the buffer that contains the vertex coordinates.
    Indices,      //!< Designates the buffer that contains the vertex indices.
    Normals,      //!< Designates the buffer that contains the normal coordinates.
    Colors,       //!< Designates the buffer that contains the color components.
    Radius,       //!< Designates the buffer that contains the radius of wire spheres.
    Centers,  //!< Designates the buffer that contains the center of c3t3 facets or the center of wire spheres, for example.
    Texture_map,        //!< Designates the buffer that contains the UV map for the texture.
    NbOfVbos      //!< Designates the size of the VBOs vector for `Edge_container`s
  };

  //!
  //! \brief The constructor.
  //! \param program is the `QOpenGLShaderProgram` that is used by this `Edge_container` `Vao`.
  //! \param indexed must be `true` if the data is indexed, `false` otherwise. If `true`, VBOs[Indices] must be filled.
  //!
  Edge_container(int program, bool indexed);

  //!
  //! \brief initGL creates the `Vbo`s and `Vao`s of this `Edge_container`.
  //! \attention It must be called within a valid OpenGL context. The `draw()` function of an item is always a safe place to call this.
  //!
  //! \todo Is it a good idea to call InitGL of each item in the scene so the developper doesn't have to worry about this in each draw() of each item ?
  //!`.
  //! \param viewer the active `Viewer_interface`.
  //!
  void initGL(Viewer_interface *viewer)  Q_DECL_OVERRIDE;

  //!
  //! \brief draw is the function that actually renders the data.
  //! \param viewer the active `Viewer_interface`.
  //! \param is_color_uniform must be `false` if `VBOs`[`Colors`] is not empty, `true` otherwise.
  //!
  void draw(CGAL::Three::Viewer_interface* viewer,
            bool is_color_uniform)  Q_DECL_OVERRIDE;

  void initializeBuffers(Viewer_interface *viewer) Q_DECL_OVERRIDE;

  /// \name Getters and Setters for the shaders parameters.
  ///
  /// Each of those depends of the `OpenGL_program_IDs` this container is using.
  /// If the shaders of this program doesn't need one, you can ignore it.
  /// The others should be filled at each `draw()` from the item.
  ///@{

  //! getter for the "plane" parameter
  QVector4D getPlane()const;
  //! getter for the "f_matrix" parameter
  QMatrix4x4 getFrameMatrix()const;
  //! getter for the "viewport" parameter
  QVector2D getViewport()const;
  //! getter for the "width" parameter
  GLfloat getWidth()const;

  //! setter for the "plane" parameter
  void setPlane(const QVector4D&);
  //! setter for the "f_matrix" parameter
  void setFrameMatrix(const QMatrix4x4&);
  //! setter for the "viewport" parameter
  void setViewport(const QVector2D&);
  //! setter for the "width" parameter
  void setWidth(const GLfloat&);
  //! setter for the "is_surface" attribute. Used in PROGRAM_C3T3_EDGES
  void setIsSurface  (const bool);
  ///@}

private:
  friend struct Edge_d;
  mutable Edge_d* d;
}; //end of class Triangle_container
}
}

#endif // EDGE_CONTAINER_H
