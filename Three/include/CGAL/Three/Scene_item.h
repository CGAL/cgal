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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Laurent RINEAU, Maxime Gimeno

#ifndef SCENE_ITEM_H
#define SCENE_ITEM_H

#include <CGAL/license/Three.h>

#include <CGAL/Three/Scene_item_config.h>
#include <CGAL/Three/Scene_interface.h>
#include <QString>
#include <QPixmap>
#include <QFont>
#include <QOpenGLBuffer>
#include <QOpenGLShader>
#include <QOpenGLVertexArrayObject>
#include <vector>
#include <CGAL/Bbox_3.h>

namespace CGAL {
namespace Three {
  class Viewer_interface;
}
}
namespace qglviewer {
  class ManipulatedFrame;
}

class QMenu;
class QKeyEvent;
namespace CGAL {
namespace Three {

class Scene_group_item;
class Viewer_interface;
//! This class represents an object in the OpenGL scene.
//! It contains all the functions called by the Scene. It
//! acts like a mix between an interface and a helper.
class SCENE_ITEM_EXPORT Scene_item : public QObject{
  Q_OBJECT
  Q_PROPERTY(QColor color READ color WRITE setColor)
  Q_PROPERTY(QString name READ name WRITE setName)
  Q_PROPERTY(bool visible READ visible WRITE setVisible)
  Q_ENUMS(RenderingMode)
  Q_PROPERTY(RenderingMode renderingMode READ renderingMode WRITE setRenderingMode)
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
  PROGRAM_CUTPLANE_SPHERES,    /** Used to render the spheres of an item with a cut plane.*/
  PROGRAM_SPHERES,             /** Used to render one or several spheres.*/
  PROGRAM_FLAT,                /** Used to render flat shading without pre computing normals*/
  PROGRAM_OLD_FLAT,            /** Used to render flat shading without pre computing normals without geometry shader*/
  NB_OF_PROGRAMS               /** Holds the number of different programs in this enum.*/
 };
  typedef CGAL::Bbox_3 Bbox;
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;
  //! \brief The default color of a scene_item.
  //!
  //! This color is the one that will be displayed if none is specified after its creation.
  static const QColor defaultColor; // defined in Scene_item.cpp

  //!\brief The Constructor.
  //!
  //! This is where the vectors of VBOs and VAOs are initialized.
  Scene_item(int buffers_size = 20, int vaos_size = 10);

  //! \brief Sets the number of isolated vertices.
  //!
  //! This number will be displayed in a warning box at loading.
  //! @see getNbIsolatedvertices
  void setNbIsolatedvertices(std::size_t nb) { nb_isolated_vertices = nb;}
  //! Getter for the number of isolated vertices.
  //! @see setNbIsolatedvertices
  std::size_t getNbIsolatedvertices() const {return nb_isolated_vertices;}
  virtual ~Scene_item();
  //! \brief Duplicates the item.
  //!
  //! Creates a new item as a copy of this one.
  virtual Scene_item* clone() const = 0;

  //! \brief Indicates if `m` is supported
  //!
  //! If it is, it will be displayed in the context menu of the item.
  virtual bool supportsRenderingMode(RenderingMode m) const = 0;
  //! Deprecated. Does nothing.
  virtual void draw() const {}
  /*! \brief The drawing function for faces.
   *
   * Draws the faces of the item in the viewer. The data
   * for the drawing is gathered in computeElements(), and is sent
   * to buffers in initializeBuffers().
   * @see computeElements()
   * @see initializeBuffers()
   */
  virtual void draw(CGAL::Three::Viewer_interface*) const  { draw(); }
  //! Deprecated. Does nothing.
  virtual void drawEdges() const { draw(); }
  /*! \brief The drawing function for the edges.
   *
   * Draws the edges and lines of the item in the viewer. The data
   * for the drawing is gathered in computeElements(), and is sent
   * to buffers in initializeBuffers().
   * @see computeElements()
   * @see initializeBuffers()
   */
  virtual void drawEdges(CGAL::Three::Viewer_interface* viewer) const { draw(viewer); }
  //! Deprecated. Does nothing.
  virtual void drawPoints() const { draw(); }
  /*! \brief The drawing function for the points.
   *
   * Draws the points of the item in the viewer. The data
   * for the drawing is gathered in computeElements(), and is sent
   * to buffers in initializeBuffers().
   * @see computeElements()
   * @see initializeBuffers()
   */
  virtual void drawPoints(CGAL::Three::Viewer_interface*) const { drawPoints(); }

  //! Called by the scene. If b is true, then this item is currently selected.
  virtual void selection_changed(bool b);

  // Functions for displaying meta-data of the item
  //!\brief Contains meta-data about the item.
  //! @returns a QString containing meta-data about the item.
  virtual QString toolTip() const = 0;
  //! \brief Contains graphical meta-data about the item.
  //! @returns a QPixmap containing graphical meta-data about the item.
  virtual QPixmap graphicalToolTip() const { return QPixmap(); }
  //! \brief Contains the font used for the data of the item.
  //! @returns a QFont containing the font used for the data of the item.
  virtual QFont font() const { return QFont(); }

  // Functions that help the Scene to compute its bbox
  //! \brief Determines if the item is finite or not.
  //!
  //! For example, a plane is not finite.
  //! If false, the BBox is not computed.
  virtual bool isFinite() const { return true; }
  //! Specifies if the item is empty or null.
  //! If true, the BBox is not computed.
  virtual bool isEmpty() const { return true; }
  //! \brief The item's bounding box.
  //!
  //! If the Bbox has never been computed, computes it and
  //! saves the result for further calls.
  //! @returns the item's bounding box.
  virtual Bbox bbox() const {
      if(!is_bbox_computed)
          compute_bbox();
      is_bbox_computed = true;
      return _bbox;
  }
  //! \brief the item's bounding box's diagonal length.
  //!
  //! If the diagonal's length has never been computed, computes it and
  //! saves the result for further calls.
  //! @returns the item's bounding box's diagonal length.
  virtual double diagonalBbox() const {
   if(!is_diag_bbox_computed)
       compute_diag_bbox();
   is_diag_bbox_computed = true;
   return _diag_bbox;
  }

  // Function about manipulation
  //! Returns true if the item has a ManipulatedFrame.
  //! @see manipulatedFrame()
  virtual bool manipulatable() const { return false; }
  //!\brief The manipulatedFrame of the item.
  //!
  //! A manipulated frame is an independant system that can be
  //! translated or rotated using the Ctrl key and the mouse.
  //!@returns the manipulatedFrame of the item.
  virtual ManipulatedFrame* manipulatedFrame() { return 0; }

  // Getters for the four basic properties
  //!Getter for the item's color.
  //! @returns the current color of the item.
  virtual QColor color() const { return color_; }
  //!Getter for the item's name.
  //! @returns the current name of the item.
  virtual QString name() const { return name_; }
  //! If the item is not visible, it is not drawn and its Bbox
  //! is ignored in the computation of the scene's.
  //! @returns the current visibility of the item.
  virtual bool visible() const { return visible_; }
  //!Getter for the item's rendering mode.
  //! @returns the current rendering mode of the item.
  //!@see RenderingMode
  virtual RenderingMode renderingMode() const { return rendering_mode; }
  //!The renderingMode's name.
  //! @returns the current rendering mode of the item as a human readable string.
  virtual QString renderingModeName() const;

  //! \brief Context menu
  //!
  //! Contains the list of the supported rendering modes,
  //! the Operations menu, actions to save or clone the item if it is supported
  //! and any contextual action for the item.
  virtual QMenu* contextMenu();

  //!Handles key press events.
  virtual bool keyPressEvent(QKeyEvent*){return false;}

  //!The group containing the item.
  //! \returns the parent group if the item is in a group
  //! \returns 0 if the item is not in a group.
  Scene_group_item* parentGroup() const;

  //!Contains the header for the table in the statistics dialog
  /*!
   * A header data is composed of 2 columns : the Categories and the titles.
   * A category is the name given to an association of titles.
   * A title is the name of a line.
   *\verbatim
   * For example,
   * Category :    | Titles| Values
   * 2 lines       |       |
   *  ____________________________
   * |             |Name   |Cube |
   * |             |_______|_____|
   * |General Info | #Edges|12   |
   * |_____________|_______|_____|
   *
   *  would be stored as follows :
   * categories = std::pair<QString,int>(QString("General Info"),2)
   * titles.append("Name");
   * titles.append("#Edges");\endverbatim
   */
  struct Header_data{
   //!Contains the name of the category of statistics and the number of lines it will contain
   QList<std::pair<QString, int> > categories;
   //!Contains the name of the lines of each category. Must be sorted from top to bottom.
   QList<QString> titles;
  };
  //!Returns a Header_data struct containing the header information.
  virtual Header_data header()const;
  //!Returns true if the item has statistics.
  virtual bool has_stats()const{return false;}
  //!Returns a QString containing the requested value for the the table in the statistics dialog
  /*! \verbatim
   * Example :
   *  ____________________________
   * |             |Name   |Cube |
   * |             |_______|_____|
   * |General Info | #Edges|12   |
   * |_____________|_______|_____|
   * compute stats(0) should return "Cube" and computeStats(1) should return QString::number(12);
   * The numbers must be coherent with the order of declaration of the titles in the header.
   * \endverbatim
   *
   */
  virtual QString computeStats(int i);

  //!Contains the number of group and subgroups containing this item.
  int has_group;

public Q_SLOTS:

  //! Notifies the program that the internal data or the properties of
  //! an item has changed, and that it must be computed again. It is
  //! important to call this function whenever the internal data is changed,
  //! or the displayed item will not be updated.
  virtual void invalidateOpenGLBuffers();
  //!Setter for the color of the item.
  virtual void setColor(QColor c) { color_ = c;}
  //!Setter for the RGB color of the item. Calls setColor(QColor).
  //!@see setColor(QColor c)
  void setRbgColor(int r, int g, int b) { setColor(QColor(r, g, b)); }
  //!Sets the name of the item.
  virtual void setName(QString n) { name_ = n; }
    //!Sets the visibility of the item.
  virtual void setVisible(bool b);
  //!Set the parent group. If `group==0`, then the item has no parent.
  //!This function is called by `Scene::changeGroup` and should not be
  //!called manually.
  virtual void moveToGroup(Scene_group_item* group);
  //!Sets the rendering mode of the item.
  //!@see RenderingMode
  virtual void setRenderingMode(RenderingMode m) { 
    if (supportsRenderingMode(m))
      rendering_mode = m; 
    Q_EMIT redraw();
  }
  //!Sets the RenderingMode to Points.
  void setPointsMode() {
    setRenderingMode(Points);
  }
  //!Sets the RenderingMode to Points.
  void setShadedPointsMode() {
    setRenderingMode(ShadedPoints);
  }
  //!Sets the RenderingMode to Wireframe.
  void setWireframeMode() {
    setRenderingMode(Wireframe);
  }

  //!Sets the RenderingMode to Flat.
  void setFlatMode() {
    setRenderingMode(Flat);
  }
  //!Set the RenderingMode to FlatPlusEdges.
  void setFlatPlusEdgesMode() {
    setRenderingMode(FlatPlusEdges);
  }
  //!Sets the RenderingMode to Gouraud.
  void setGouraudMode() {
    setRenderingMode(Gouraud);
  }
  //!Sets the RenderingMode to PointsPlusNormals.
  void setPointsPlusNormalsMode(){
    setRenderingMode(PointsPlusNormals);
  }
  
  //!Emits an aboutToBeDestroyed() signal.
  //!Override this function to delete what needs to be deleted on destruction.
  //!This might be needed as items are not always deleted right away by Qt and this behaviour may cause a simily
  //!memory leak, for example when multiple items are created at the same time.
  virtual void itemAboutToBeDestroyed(Scene_item*);

  //!Selects a point through raycasting.
  virtual void select(double orig_x,
                      double orig_y,
                      double orig_z,
                      double dir_x,
                      double dir_y,
                      double dir_z);



Q_SIGNALS:
  //! Is emitted to notify a change in the item's data.
  void itemChanged();
  //! Is emitted when the item is shown to notify a change in the item's visibility.
  //! Typically used to update the scene's bbox;
  void itemVisibilityChanged();
  //! Is emitted to notify that the item is about to be deleted.
  void aboutToBeDestroyed();
  //! Is emitted to require a new display.
  void redraw();

protected:
  //!Holds the BBox of the item
  mutable Bbox _bbox;
  mutable double _diag_bbox;
  mutable bool is_bbox_computed;
  mutable bool is_diag_bbox_computed;
  virtual void compute_bbox()const{}
  virtual void compute_diag_bbox()const;
  // The four basic properties
  //!The name of the item.
  QString name_;
  //!The color of the item.
  QColor color_;
  //!The visibility of the item.
  bool visible_;
  //!The parent group, or 0 if the item is not in a group.
  Scene_group_item* parent_group;
  //!Specifies if the item is currently selected.
  bool is_selected;
  //! Holds the number of vertices that are not linked to the polyhedron from the OFF
  //! file.
  std::size_t nb_isolated_vertices;
  /*! Decides if the draw function must call initializeBuffers() or not. It is set
   * to true in the end of initializeBuffers() and to false in invalidateOpenGLBuffers(). The need of
   * this boolean comes from the need of a context from the OpenGLFunctions used in
   * initializeBuffers().
   * @see initializeBuffers()
   * @see invalidateOpenGLBuffers()
   */
  mutable bool are_buffers_filled;
  //!The rendering mode of the item.
  //!@see RenderingMode
  RenderingMode rendering_mode;
  //!The default context menu.
  QMenu* defaultContextMenu;
  /*! Contains the previous RenderingMode.
   * This is used to determine if invalidateOpenGLBuffers should be called or not
   * in certain cases.
   * @see invalidateOpenGLBuffers()*/
  RenderingMode prev_shading;
  /*! \todo replace it by RenderingMode().
   * \brief Contains the current RenderingMode.
   *
   * This is used to determine if invalidateOpenGLBuffers should be called or not
   * in certain cases.
   * @see invalidateOpenGLBuffers()*/
  RenderingMode cur_shading;
  //!Contains the size of the vector of VBOs
  int buffersSize;
  //!Contains the size of the map of VAOs
  int vaosSize;
  //!Contains the VBOs
  mutable std::vector<QOpenGLBuffer> buffers;
  /*! Contains the VAOs.
   */
  std::vector<QOpenGLVertexArrayObject*> vaos;
  //!Adds a VAO to the Map.
  void addVaos(int i)
  {
      QOpenGLVertexArrayObject* n_vao = new QOpenGLVertexArrayObject();
      vaos[i] = n_vao;
  }

  /*! Fills the VBOs with data. Must be called after each call to #computeElements().
   * @see compute_elements()
   */
  void initializeBuffers(){}

  /*! Collects all the data for the shaders. Must be called in #invalidateOpenGLBuffers().
   * @see invalidateOpenGLBuffers().
   */
  void computeElements(){}
  /*! Passes all the uniform data to the shaders.
   * According to program_name, this data may change.
   */
  void attribBuffers(CGAL::Three::Viewer_interface*, int program_name) const;

  /*! Compatibility function. Calls `viewer->getShaderProgram()`. */
  virtual QOpenGLShaderProgram* getShaderProgram(int name , CGAL::Three::Viewer_interface *viewer = 0) const;
}; // end class Scene_item
}
}

#include <QMetaType>
Q_DECLARE_METATYPE(CGAL::Three::Scene_item*)

#endif // SCENE_ITEM_H
