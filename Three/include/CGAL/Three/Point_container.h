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

#ifndef POINT_CONTAINER_H
#define POINT_CONTAINER_H

#include <CGAL/license/Three.h>


#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Primitive_container.h>


#ifdef demo_framework_EXPORTS
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_EXPORT
#else
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_IMPORT
#endif
struct Point_d;
namespace CGAL {
namespace Three {

//!
//! \brief The Point_container struct wraps the OpenGL data for drawing lines.
//!
struct DEMO_FRAMEWORK_EXPORT Point_container :public Primitive_container
{

  //!
  //! \brief The vbosName enum
  //!
  //! Holds the `Vbo` Ids of this container.
  //!
  enum vbosName {
    Vertices = 0, //!< Designates the buffer that contains the vertex coordinates.
    Indices,      //!< Designates the buffer that contains the vertex indices.
    Colors,       //!< Designates the buffer that contains the color components.
    NbOfVbos      //!< Designates the size of the VBOs vector for `Point_container`s
  };

  //!
  //! \brief The constructor.
  //! \param program is the `QOpenGLShaderProgram` that is used by this `Point_container` `Vao`.
  //! \param indexed must be `true` if the data is indexed, `false` otherwise. If `true`, VBOs[Indices] must be filled.
  //!
  Point_container(int program, bool indexed);

  //!
  //! \brief initGL creates the `Vbo`s and `Vao`s of this `Point_container`.
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
//! setter for the "f_matrix" parameter
  void setFrameMatrix(const QMatrix4x4&);
  ///@}

private:
  friend struct ::Point_d;
  mutable ::Point_d* d;
}; //end of class Triangle_container
}
}

#endif // POINT_CONTAINER_H
