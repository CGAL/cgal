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

#ifndef TRIANGLE_CONTAINER_H
#define TRIANGLE_CONTAINER_H

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
struct Tri_d;
namespace CGAL {
namespace Three {

//! \brief The Triangle_container struct wraps the OpenGL data for drawing triangles.
struct DEMO_FRAMEWORK_EXPORT Triangle_container :public Primitive_container
{

  //! \brief The vbosName enum
  //!
  //! Holds the `Vbo` Ids of this container.
  //!
  enum vbosName {
    Flat_vertices = 0,  //!< Designates the buffer that contains the flat vertex coordinates (not indexed).
    Smooth_vertices,    //!< Designates the buffer that contains the smooth vertex coordinates (indexed).
    Vertex_indices,     //!< Designates the buffer that contains the indices for the smooth vertices.
    Flat_normals,       //!< Designates the buffer that contains the normals for the flat vertices.
    Smooth_normals,     //!< Designates the buffer that contains the normals for the smooth vertices.
    Facet_barycenters,  //!< Designates the buffer that contains the barycenters of the c3t3 facets or the center of the spheres.
    Radius,             //!< Designates the buffer that contains the radius of the spheres.
    VColors,            //!< Designates the buffer that contains the colors of the smooth vertices.
    FColors,            //!< Designates the buffer that contains the colors of the flat vertices.
    NbOfVbos            //!< Designates the size of the VBOs vector for `Triangle_container`s
  };

  //!
  //! \brief The constructor.
  //! \param program is the `QOpenGLShaderProgram` that is used by this `Triangle_container` `Vao`.
  //! \param indexed must be `true` if the data is indexed, `false` otherwise. If `true`, `VBOs`[`Vertex_indices`] must be filled.
  //!
  Triangle_container(int program, bool indexed);

  //!
  //! \brief initGL creates the Vbos and Vaos of this `Triangle_container`.
  //! \attention It must be called within a valid OpenGL context. The `draw()` function of an item is always a safe place to call this.
  //!
  //! \param viewer the active `Viewer_interface`.
  //!
  void initGL(CGAL::Three::Viewer_interface* viewer) Q_DECL_OVERRIDE;

  //!
  //! \brief draw is the function that actually renders the data.
  //! \param viewer the active `Viewer_interface`.
  //! \param is_color_uniform must be `false` if the color buffers are not empty, `true` otherwise.
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
  
  //! getter for the "shrink_factor" parameter
  float getShrinkFactor();
  //! getter for the "plane" parameter
  QVector4D getPlane();
  //! getter for the "alpha" parameter
  float getAlpha();
//! setter for the "shrink_factor" parameter
  void setShrinkFactor(const float&);
  //! setter for the "plane" parameter
  void setPlane       (const QVector4D&);
  //! setter for the "alpha" parameter
  void setAlpha       (const float&);
  ///@}

  //drawing variables


private:
  friend struct Tri_d;
  mutable Tri_d* d;

}; //end of class Triangle_container

}
}

#endif // TRIANGLE_CONTAINER_H
