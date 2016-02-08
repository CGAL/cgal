// Copyright (c) 2012-2015  GeometryFactory Sarl (France)
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
// Author(s)     : Laurent RINEAU, Maxime Gimeno

#ifndef VIEWER_INTERFACE_H
#define VIEWER_INTERFACE_H

#include <QGLViewer/qglviewer.h>
#include <QWidget>
#include <QPoint>
#include <QOpenGLFunctions_2_1>
#include <CGAL/Qt/CreateOpenGLContext.h>
// forward declarations
class QWidget;
class QMouseEvent;
class QKeyEvent;
class QOpenGLShaderProgram;

//! \file Viewer_interface.h
#include <CGAL/Three/Viewer_config.h> // for VIEWER_EXPORT
namespace CGAL{
namespace Three{
class Scene_draw_interface;
//! Base class to interact with the viewer from the plugins, the items and the scene.
class VIEWER_EXPORT Viewer_interface : public QGLViewer, public QOpenGLFunctions_2_1 {

  Q_OBJECT

public:
 /*!
   * \brief The OpenGL_program_IDs enum
   * This enum holds the OpenGL programs IDs that are given to getShaderProgram() and attrib_buffers().
   *@see getShaderProgram
   * @see attrib_buffers
   */
  enum OpenGL_program_IDs
  {
   PROGRAM_WITH_LIGHT = 0,      /** Used to render a surface or edge affected by the light. It uses a per fragment lighting model, and renders brighter the selected item.*/
   PROGRAM_WITHOUT_LIGHT,       /** Used to render a polygon edge or points. It renders in a uniform color and is not affected by light. It renders the selected item in black.*/
   PROGRAM_NO_SELECTION,        /** Used to render a polyline or a surface that is not affected by light, like a cutting plane. It renders in a uniform color that does not change with selection.*/
   PROGRAM_WITH_TEXTURE,        /** Used to render a textured polyhedron. Affected by light.*/
   PROGRAM_PLANE_TWO_FACES,     /** Used to render a two-faced plane. The two faces have a different color. Not affected by light.*/
   PROGRAM_WITH_TEXTURED_EDGES, /** Used to render the edges of a textured polyhedorn. Not affected by light.*/
   PROGRAM_INSTANCED,           /** Used to display instanced rendered spheres.Affected by light.*/
   PROGRAM_INSTANCED_WIRE,      /** Used to display instanced rendered wired spheres. Not affected by light.*/
   PROGRAM_C3T3,                /** Used to render a c3t3_item. It discards any fragment on a side of a plane, meaning that nothing is displayed on this side of the plane. Affected by light.*/
   PROGRAM_C3T3_EDGES,          /** Used to render the edges of a c3t3_item. It discards any fragment on a side of a plane, meaning that nothing is displayed on this side of the plane. Not affected by light.*/
   NB_OF_PROGRAMS               /** Holds the number of different programs in this enum.*/
  };

  Viewer_interface(QWidget* parent) : QGLViewer(CGAL::Qt::createOpenGLContext(), parent) {}
  virtual ~Viewer_interface() {}

  //! Sets the scene for the viewer.
  virtual void setScene(CGAL::Three::Scene_draw_interface* scene) = 0;
  //! @returns the antialiasing state.
  virtual bool antiAliasing() const = 0;
  
  // Those two functions are defined in Viewer.cpp
  //!Sets the position and orientation of a frame using a QString.
  //!@returns true if it worked.
  static bool readFrame(QString, qglviewer::Frame&);
  //!@returns a QString contining the position and orientation of a frame.
  static QString dumpFrame(const qglviewer::Frame&);
  //! @returns the fastDrawing state.
  virtual bool inFastDrawing() const = 0;

  /*! Passes all the uniform data to the shaders.
   * According to program_name, this data may change.
   * @see OpenGL_program_IDs
   */
  virtual void attrib_buffers(int program_name) const = 0;

  /*! Returns a program according to name.
   * If the program does not exist yet, it is created and stored in shader_programs.
   * @see OpenGL_program_IDs
   * @returns a pointer to the corresponding program.*/
  virtual QOpenGLShaderProgram* getShaderProgram(int name) const = 0;

  //!Allows OpenGL 2.1 context to get access to glDrawArraysInstanced.
  typedef void (APIENTRYP PFNGLDRAWARRAYSINSTANCEDARBPROC) (GLenum mode, GLint first, GLsizei count, GLsizei primcount);
  //!Allows OpenGL 2.1 context to get access to glVertexAttribDivisor.
  typedef void (APIENTRYP PFNGLVERTEXATTRIBDIVISORARBPROC) (GLuint index, GLuint divisor);
  //!Allows OpenGL 2.1 context to get access to glVertexAttribDivisor.
  typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTURE2DEXTPROC) (GLuint target, GLuint attachment, GLuint textarget, GLuint texture, GLint level);
  //!Allows OpenGL 2.1 context to get access to gkFrameBufferTexture2D.
  PFNGLDRAWARRAYSINSTANCEDARBPROC glDrawArraysInstanced;
  //!Allows OpenGL 2.1 context to get access to glVertexAttribDivisor.
  PFNGLVERTEXATTRIBDIVISORARBPROC glVertexAttribDivisor;
  //!Allows OpenGL 2.1 context to get access to gkFrameBufferTexture2D.
  PFNGLFRAMEBUFFERTEXTURE2DEXTPROC glFramebufferTexture2D;
  //!@returns true if glVertexAttribDivisor, and glDrawArraysInstanced are found.
  //! Used by the items to avoid SEGFAULT.
  bool extension_is_found;
  //!The matrix used for the picking.
  GLfloat pickMatrix_[16];
  //!Sets the binding for SHIFT+LEFT CLICK to SELECT (initially used in Scene_polyhedron_selection_item.h)
  void setBindingSelect()
  {
#if QGLVIEWER_VERSION >= 0x020501
    setMouseBinding(::Qt::ShiftModifier, ::Qt::LeftButton, SELECT);
#else
    setMouseBinding(::Qt::SHIFT + ::Qt::LeftButton, SELECT);
#endif
  }
  //!Sets the binding for SHIFT+LEFT CLICK to NO_CLICK_ACTION (initially used in Scene_polyhedron_selection_item.h)
  void setNoBinding()
  {
#if QGLVIEWER_VERSION >= 0x020501
    setMouseBinding(::Qt::ShiftModifier, ::Qt::LeftButton, NO_CLICK_ACTION);
#else
    setMouseBinding(::Qt::SHIFT + ::Qt::LeftButton, NO_CLICK_ACTION);
#endif
  }

Q_SIGNALS:
  //!Defined automatically in moc.
  void selected(int);
  //!Defined automatically in moc.
  void requestContextMenu(QPoint global_pos);
  //!Defined automatically in moc.
  void selectedPoint(double, double, double);
  //!Defined automatically in moc.
  void selectionRay(double, double, double, double, double, double);

public Q_SLOTS:
  //! Sets the antialiasing to true or false.
  virtual void setAntiAliasing(bool b) = 0;
  //! If b is true, facets will be ligted from both internal and external sides.
  //! If b is false, only the side that is exposed to the light source will be lighted.
  virtual void setTwoSides(bool b) = 0;
  //! If b is true, some items are displayed in a simplified version when moving the camera.
  //! If b is false, items display is never altered, even when moving.
  virtual void setFastDrawing(bool b) = 0;
  //! Make the camera turn around.
  virtual void turnCameraBy180Degres() = 0;
  //! @returns a QString containing the position and orientation of the camera.
  virtual QString dumpCameraCoordinates() = 0;
//!Moves the camera to the new coordinates (position and orientation) through an animation.
  virtual bool moveCameraToCoordinates(QString, 
                                       float animation_duration = 0.5f) = 0;

}; // end class Viewer_interface
}
}
#endif // VIEWER_INTERFACE_H
