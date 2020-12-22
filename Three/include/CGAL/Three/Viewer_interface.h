// Copyright (c) 2012-2015  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU, Maxime Gimeno

#ifndef VIEWER_INTERFACE_H
#define VIEWER_INTERFACE_H

#include <CGAL/license/Three.h>

#include <QMap>
#include <CGAL/Qt/qglviewer.h>
#include <QWidget>
#include <QPoint>
#include <QOpenGLFunctions>
#include <QOpenGLFunctions_4_3_Core>
#include <CGAL/Qt/CreateOpenGLContext.h>
// forward declarations
class QWidget;
class QImage;
class QMouseEvent;
class QKeyEvent;
class QOpenGLShaderProgram;
class QOpenGLFramebufferObject;
class TextRenderer;
class TextListItem;

//! \file Viewer_interface.h
#include <CGAL/Three/Viewer_config.h> // for VIEWER_EXPORT
namespace CGAL{
namespace Three{
class Scene_draw_interface;
class Scene_item;
//! Base class to interact with the viewer from the plugins, the items and the scene.
class VIEWER_EXPORT Viewer_interface : public CGAL::QGLViewer{

  Q_OBJECT

public:
  /*!
    * \brief The OpenGL_program_IDs enum
    *
    * This enum holds the OpenGL programs IDs that are given to getShaderProgram() and attribBuffers().
    * @see getShaderProgram
    * @see attribBuffers
    */
  enum OpenGL_program_IDs
  {
   PROGRAM_WITH_LIGHT = 0,      //! Used to render a surface or an edge affected by the light. It uses a per fragment lighting model, and renders the selected item brighter.
   PROGRAM_WITHOUT_LIGHT,       //! Used to render a polyhedron edge or points. It renders in a uniform color and is not affected by light. \attention It renders the selected item in black.
   PROGRAM_NO_SELECTION,        //! Used to render a polyline or a surface that is not affected by light, like a cutting plane. It renders in a uniform color that does not change with selection.
   PROGRAM_WITH_TEXTURE,        //! Used to render a textured polyhedron. Affected by light.
   PROGRAM_PLANE_TWO_FACES,     //! Used to render a two-faced plane. The two faces have a different color. Not affected by light.
   PROGRAM_WITH_TEXTURED_EDGES, //! Used to render the edges of a textured polyhedron. Not affected by light.
   PROGRAM_INSTANCED,           //! Used to display instanced rendered spheres.Affected by light.
   PROGRAM_INSTANCED_WIRE,      //! Used to display instanced rendered wired spheres. Not affected by light.
   PROGRAM_C3T3,                //! Used to render a c3t3_item. It discards any fragment on a side of a plane, meaning that nothing is displayed on this side of the plane. Affected by light.
   PROGRAM_C3T3_EDGES,          //! Used to render the edges of a c3t3_item. It discards any fragment on a side of a plane, meaning that nothing is displayed on this side of the plane. Not affected by light.
   PROGRAM_CUTPLANE_SPHERES,    //! Used to render the spheres of an item with a cut plane.
   PROGRAM_SPHERES,             //! Used to render one or several spheres.
   PROGRAM_DARK_SPHERES,        //! Used to render one or several spheres without light (for picking for example).
   PROGRAM_FLAT,                /** Used to render flat shading without pre computing normals*/
   PROGRAM_OLD_FLAT,            /** Used to render flat shading without pre computing normals without geometry shader*/
   PROGRAM_SOLID_WIREFRAME,     //! Used to render edges with width superior to 1.
   PROGRAM_NO_INTERPOLATION,   //! Used to render faces without interpolating their color.
   PROGRAM_HEAT_INTENSITY,      //! Used to render special item in Display_property_plugin
   NB_OF_PROGRAMS               //! Holds the number of different programs in this enum.
  };

 //! \brief The viewer's QPainter
 //!
 //! The painter is the element that draws everything on screen,
 //! but you should only need this if you want to draw 2D things
 //! on top of the scene, like a selection rectangle.
 //! See <a href="https://doc.qt.io/qt-5/qpainter.html">QPainter's documentation </a> for details.
 virtual QPainter *getPainter() =0;



  //! \brief Tests if an id should be displayed or not.
  //! \param x, y, z the coordinates of the id's position.
  //! \return true if the ID should be displayed.
  virtual bool testDisplayId(double x, double y, double z) = 0;
  //! \brief Updates the item's displayed ids.
  //!
  //! Call this after the mesh or its ids have changed.
  virtual void updateIds(CGAL::Three::Scene_item *) = 0;
  //! \brief Specifies if the items ids are being displayed.
  //!
  //! \returns true if the primitive ids are currently displayed
  virtual bool hasText() const { return false; }
  //! \brief Constructor
  //!
  //! Creates a valid context for OpenGL ES 2.0.
  //! \param parent the parent widget. It usually is the MainWindow.
  Viewer_interface(QWidget* parent) : CGAL::QGLViewer(parent) {}
  //!
  //! \brief Constructor for the secondary viewers.
  //!
  //! \param parent the parent widget. It usually is the MainWindow.
  //! \param shared_widget the main viewer of the Application. This will share the
  //!  context and allow synchronized rendering of multiple views.
  //!
  Viewer_interface(QWidget* parent, QOpenGLWidget* shared_widget)
    : QGLViewer(shared_widget->context(),parent){}
  virtual ~Viewer_interface() {}

  //! \brief Sets the scene for the viewer.
  virtual void setScene(CGAL::Three::Scene_draw_interface* scene) = 0;
  //! \brief The antialiasing state.
  //!
  //! @returns true if the antialiasing is activated.
  virtual bool antiAliasing() const = 0;

  // Those two functions are defined in Viewer.cpp
  //! \brief Sets the position and orientation of a frame using a QString.
  //! \param s is usually gotten by dumpFrame() and is of the form "Px Py Pz O1 O2 O3 O4 ", with
  //! - Px to Py : the new position coordinates,
  //! - O1 to O3 : axis coordinate *sin(angle/2)
  //! - O4 cos(angle/2).
  //! \param frame is the frame that will be moved
  //! @returns true if it worked.
  //! @see moveCameraToCoordinates()
  static bool readFrame(QString s, CGAL::qglviewer::Frame& frame);
  //! \brief Gives information about a frame.
  //! @see readFrame
  //! @see dumpCameraCoordinates()
  //!@returns a QString containing the position and orientation of a frame.
  static QString dumpFrame(const CGAL::qglviewer::Frame&);
  //! \brief The fastDrawing state.
  //!
  //! In fast drawing mode, some items will be simplified while the camera is moving
  //! to gain in performance.
  //! @returns the fastDrawing state.
  virtual bool inFastDrawing() const = 0;
  //! \brief The drawWithNames state.
  //!
  //! In draw with name mode, the scene is not displayed, but a
  //! \a name is given to each Scene_item. It is used for picking.
  //! @returns true if the viewer is drawing with names.
  virtual bool inDrawWithNames() const = 0;

  //! \brief Passes all the uniform data to the shaders.
  //!
  //! According to program_name, this data may change.
  //! This should be called in every Scene_item::draw() call.
  //! @see OpenGL_program_IDs
  //!
  virtual void attribBuffers(int program_name) const = 0;
  /*! Enables the clipping box. Each Vector4 of `box` contains the equation of a plane of the clipping box.
   * Everything that is located on the positive side of one of those planes will not be displayed.
   * @see disableCLippingBox()
   */
  virtual void enableClippingBox(QVector4D box[6])=0;

  /*!
   * Disables the clipping box. The six clipping planes will be ignored.
   * @see enableClippingBox()
   */
  virtual void disableClippingBox()= 0;

  //! \brief Returns a program according to name.
  //!
  //! If the program does not exist yet, it is created and stored in shader_programs.
  //! @see OpenGL_program_IDs
  //! @returns a pointer to the corresponding program.
  virtual QOpenGLShaderProgram* getShaderProgram(int name) const = 0;

  //!\brief TextRenderer is used to display text on the screen.
  //!
  //! The textRenderer uses the painter to display 2D text over the 3D Scene.
  //! \returns the viewer's TextRender
  virtual TextRenderer* textRenderer() = 0;
  //!Allows OpenGL ES 2.0 context to get access to glDrawArraysInstanced.
  typedef void (APIENTRYP PFNGLDRAWARRAYSINSTANCEDARBPROC) (GLenum mode, GLint first, GLsizei count, GLsizei primcount);
  //!Allows OpenGL ES 2.0 context to get access to glVertexAttribDivisor.
  typedef void (APIENTRYP PFNGLVERTEXATTRIBDIVISORARBPROC) (GLuint index, GLuint divisor);
  //!Allows OpenGL ES 2.0 context to get access to glVertexAttribDivisor.
  typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTURE2DEXTPROC) (GLuint target, GLuint attachment, GLuint textarget, GLuint texture, GLint level);

  PFNGLDRAWARRAYSINSTANCEDARBPROC glDrawArraysInstanced;
  PFNGLVERTEXATTRIBDIVISORARBPROC glVertexAttribDivisor;
  PFNGLFRAMEBUFFERTEXTURE2DEXTPROC glFramebufferTexture2D;

  //! \brief Used by the items to avoid SEGFAULT.
  //!@returns true if glVertexAttribDivisor, and glDrawArraysInstanced are found.
  virtual bool isExtensionFound() = 0;
  //!\brief Allows to perform picking from the keyboard and mouse
  //!
  //! Sets the combination SHIFT+LEFT CLICK to perform a selection on the screen.
  //! This is used to perform picking.
  virtual void setBindingSelect() = 0;
  //!\brief Disable the picking from the keyboard and mouse
  //!
  //! Unbinds the combination SHIFT+LEFT CLICK. It allows to
  //! avoid conflicts in the selection_tool, for example.
  virtual void setNoBinding() = 0 ;

  //!
  //! If this mode is ON, the viewer will display the content of `staticImage()` instead
  //! of drawing the cene. This is used when drawing 2D lines over the viewer.
  //! @see `staticImage()`
  //! @see `setStaticImage()`
  virtual void set2DSelectionMode(bool) = 0;

  //!
  //! Setter for the image to be displayed in 2D selection mode.
  //!
  virtual void setStaticImage(QImage image)=0;

  //! Returns the static image to be displayed in 2D selection mode.
  virtual const QImage& staticImage() const = 0;

  //!The number of passes that are performed for the scene transparency.
  //! Customizable from the MainWindow or the SubViewer menu.
  virtual float total_pass() = 0;
Q_SIGNALS:
  //!Emit this to signal that the `id`th item has been picked.
  void selected(int id);
  //!Emit this to require a contextual menu to appear at `global_pos`.
  void requestContextMenu(QPoint global_pos);
  //!Emit this to signal that the point at (`x`, `y`, `z`) has been picked.
  void selectedPoint(double x, double y, double z);
  //!Emit this to request the currently selected item to perform a selection based on an AABB_Tree and a raycasting.
  void selectionRay(double sx, double sy, double sz, double tx, double ty, double tz);

public Q_SLOTS:
  //! Sets the antialiasing to true or false.
  //! @see antiAliasing()
  virtual void setAntiAliasing(bool b) = 0;
  //! If b is true, faces will be ligted from both internal and external side.
  //! If b is false, only the side that is exposed to the light source will be lighted.
  virtual void setTwoSides(bool b) = 0;
  //! If b is true, then a special color mask is applied to points and meshes to differenciate
  //! front-faced and back-faced elements.
  virtual void setBackFrontShading(bool b) =0;
  //! \brief Sets the fast drawing mode
  //! @see inFastDrawing()
  virtual void setFastDrawing(bool b) = 0;
  //! Makes the camera turn around.
  virtual void turnCameraBy180Degres() = 0;
  //! @returns a QString containing the position and orientation of the camera.
  //! @see dumpFrame()
  virtual QString dumpCameraCoordinates() = 0;
//! \brief Moves the camera to the new coordinates.
//!
//! The movement is performed through an animation.
//! \param target is usually gotten by dumpCameraCoordinates() and is of the form "Px Py Pz O1 O2 O3 O4 ", with
//! - Px to Py : the new position coordinates,
//! - O1 to O3 : axis coordinate *sin(angle/2)
//! - O4 cos(angle/2).
//! \param animation_duration is the duration of the animation of the movement.
  virtual bool moveCameraToCoordinates(QString target,
                                       float animation_duration = 0.5f) = 0;
  //!
  //! Setter for the orthogonal projection of the viewer.
  //!
  virtual void SetOrthoProjection( bool b) =0;
public:

  //! Gives acces to recent openGL(4.3) features, allowing use of things like
  //! Geometry Shaders or Depth Textures.
  //! @returns a pointer to an initialized  QOpenGLFunctions_4_3_Core if `isOpenGL_4_3()` is `true`
  //! @returns nullptr if `isOpenGL_4_3()` is `false`
  virtual QOpenGLFunctions_4_3_Core* openGL_4_3_functions() = 0;
  //! getter for point size under old openGL context;
  virtual const GLfloat& getGlPointSize()const = 0;
  //! setter for point size under old openGL context;
  virtual void setGlPointSize(const GLfloat& p) = 0;
  virtual void setCurrentPass(int pass) = 0;
  virtual void setDepthWriting(bool writing_depth) = 0;
  virtual void setDepthPeelingFbo(QOpenGLFramebufferObject* fbo) = 0;

  virtual int currentPass()const = 0;
  virtual bool isDepthWriting()const = 0;
  virtual QOpenGLFramebufferObject* depthPeelingFbo() = 0;
  virtual void makeCurrent() = 0;
  virtual QVector4D* clipBox() const =0;
  virtual bool isClipping() const = 0;
  //!  A vector indicating the scaling factors to apply to the scene when displaying it.
  //!  It can be useful when a scene is very large along one of it's coordinates, making it hard to visualize it.
  virtual const QVector3D& scaler() const = 0;
}; // end class Viewer_interface
}
}
#endif // VIEWER_INTERFACE_H
