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

#ifndef SCENE_ITEM_H
#define SCENE_ITEM_H
#include <CGAL/Three/Scene_item_config.h>
#include <CGAL/Three/Scene_interface.h>
#include <QString>
#include <QPixmap>
#include <QFont>
#include <QOpenGLBuffer>
#include <QOpenGLShader>
#include <QOpenGLVertexArrayObject>
#include <vector>
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

class Viewer_interface;
//! This class represents an object in the OpenGL scene
class SCENE_ITEM_EXPORT Scene_item : public QObject {
  Q_OBJECT
  Q_PROPERTY(QColor color READ color WRITE setColor)
  Q_PROPERTY(QString name READ name WRITE setName)
  Q_PROPERTY(bool visible READ visible WRITE setVisible)
  Q_ENUMS(RenderingMode)
  Q_PROPERTY(RenderingMode renderingMode READ renderingMode WRITE setRenderingMode)
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
  typedef CGAL::Three::Scene_interface::Bbox Bbox;
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;
  //! The default color of a scene_item.
  static const QColor defaultColor; // defined in Scene_item.cpp

  //!The default Constructor.
  /*!
   * Initializes the number of VBOs to 20 and VAOs to 10 and creates them.
   */
  Scene_item()
    : name_("unamed"),
      color_(defaultColor),
      visible_(true),
      are_buffers_filled(false),
      rendering_mode(FlatPlusEdges),
      defaultContextMenu(0),
      buffersSize(20),
      vaosSize(10),
      vaos(10)
  {
      is_bbox_computed = false;
      is_monochrome = true;
      for(int i=0; i<vaosSize; i++)
      {
          addVaos(i);
          vaos[i]->create();
      }

      buffers.reserve(buffersSize);
      for(int i=0; i<buffersSize; i++)
      {
          QOpenGLBuffer n_buf;
          buffers.push_back(n_buf);
          buffers[i].create();
      }
      nb_isolated_vertices = 0;
      has_group = 0;
  }
  //!The Constructor.
  /*!
   * Initializes the number of VBOs and VAOs and creates them.
   */
  Scene_item(int buffers_size, int vaos_size)
    : name_("unamed"),
      color_(defaultColor),
      visible_(true),
      are_buffers_filled(false),
      rendering_mode(FlatPlusEdges),
      defaultContextMenu(0),
      buffersSize(buffers_size),
      vaosSize(vaos_size),
      vaos(vaos_size)
  {
      is_bbox_computed = false;
      is_monochrome = true;
      for(int i=0; i<vaosSize; i++)
      {
          addVaos(i);
          vaos[i]->create();
      }

      buffers.reserve(buffersSize);
      for(int i=0; i<buffersSize; i++)
      {
          QOpenGLBuffer n_buf;
          buffers.push_back(n_buf);
          buffers[i].create();
      }
      nb_isolated_vertices = 0;
      has_group = 0;
  }
  //! Setter for the number of isolated vertices.
  void setNbIsolatedvertices(std::size_t nb) { nb_isolated_vertices = nb;}
  //! Getter for the number of isolated vertices.
  std::size_t getNbIsolatedvertices() const {return nb_isolated_vertices;}
  //!The destructor. It is virtual as the item is virtual.
  virtual ~Scene_item();
  //! Creates a new item as a copy of the sender. Must be overloaded.
  virtual Scene_item* clone() const = 0;

  //! Indicates if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const = 0;
  //! Deprecated. Does nothing.
  virtual void draw() const {}
  /*! \brief The drawing function.
   * Draws the facets of the item in the viewer using OpenGL functions. The data
   * for the drawing is gathered in compute_elements(), and is sent
   * to buffers in initialize_buffers().
   * @see compute_elements()
   * @see initialize_buffers()
   */
  virtual void draw(CGAL::Three::Viewer_interface*) const  { draw(); }
  //! Deprecated. Does nothing.
  virtual void draw_edges() const { draw(); }
  /*! \brief The drawing function.
   * Draws the edges and lines of the item in the viewer using OpenGL functions. The data
   * for the drawing is gathered in compute_elements(), and is sent
   * to buffers in initialize_buffers().
   * @see compute_elements()
   * @see initialize_buffers()
   */
  virtual void draw_edges(CGAL::Three::Viewer_interface* viewer) const { draw(viewer); }
  //! Deprecated. Does nothing.
  virtual void draw_points() const { draw(); }
  /*! \brief The drawing function.
   * Draws the points of the item in the viewer using OpenGL functions. The data
   * for the drawing is gathered in compute_elements(), and is sent
   * to buffers in initialize_buffers().
   * @see compute_elements()
   * @see initialize_buffers()
   */
  virtual void draw_points(CGAL::Three::Viewer_interface*) const { draw_points(); }

  //! Draws the splats of the item in the viewer using GLSplat functions.
  virtual void draw_splats() const {}
  //! Draws the splats of the item in the viewer using GLSplat functions.
  virtual void draw_splats(CGAL::Three::Viewer_interface*) const {draw_splats();}

  //! Specifies which data must be updated when selection has changed.
  //! Must be overloaded.
  virtual void selection_changed(bool);

  // Functions for displaying meta-data of the item
  //! @returns a QString containing meta-data about the item.
  //! Must be overloaded.
  //! Data is :Number of vertices, Number of edges, Number of facets,
  //! Volume, Area, Bounding box limits and Number of isolated points.
  virtual QString toolTip() const = 0;
  //! @returns a QPixmap containing graphical meta-data about the item.
  virtual QPixmap graphicalToolTip() const { return QPixmap(); }
  //! @returns a QFont containing the font used for the data of the item.
  virtual QFont font() const { return QFont(); }

  // Functions that help the Scene to compute its bbox
  //! If isFinite() returns false, the BBox is not computed.
  virtual bool isFinite() const { return true; }
  //! Specifies if the item is empty or null.
  virtual bool isEmpty() const { return true; }
  //!@returns the item's bounding box.
  virtual Bbox bbox() const {
      if(!is_bbox_computed)
          compute_bbox();
      is_bbox_computed = true;
      return _bbox;
  }
  // Function about manipulation
  //! Decides if the item can have a ManipulatedFrame.
  virtual bool manipulatable() const { return false; }
  //!@returns the manipulatedFrame of the item.
  virtual ManipulatedFrame* manipulatedFrame() { return 0; }

  // Getters for the four basic properties
  //! @returns the current color of the item.
  virtual QColor color() const { return color_; }
  //! @returns the current name of the item.
  virtual QString name() const { return name_; }
  //! @returns the current visibility of the item.
  virtual bool visible() const { return visible_; }
  //! @returns the current rendering mode of the item.
  //!@see RenderingMode
  virtual RenderingMode renderingMode() const { return rendering_mode; }
  //! @returns the current rendering mode of the item as a human readable string.
  virtual QString renderingModeName() const;

  //! Context menu
  virtual QMenu* contextMenu();

  //!Handles key press events.
  virtual bool keyPressEvent(QKeyEvent*){return false;}

  //!Contains the header for the table in the statistics dialog
  /*!
   * A header data is composed of 2 columns : the Categories and the titles.
   * A category is the name given to an association of titles.
   * A title is the name of a line.
   *
   * For example,
   * Category :    | Titles| Values
   * 2 lines       |       |
   *  ____________________________
   * |             |Name   |Cube |
   * |             |_______|_____|
   * |General Info |\#Edges|12   |
   * |_____________|_______|_____|
   *
   *  would be stored as follows :
   * categories = std::pair<QString,int>(QString("General Info"),2)
   * titles.append("Name");
   * titles.append("#Edges");
   */
  struct Header_data{
   //!Contains the name of the category of statistics and the number of lines it will contain
   QList<std::pair<QString, int> > categories;
   //!Contains the name of the lines of each category. Must be sorted from top to bottom.
   QList<QString> titles;
  };
  //!Returns a Header_data struct containing the header information.
  virtual Header_data header()const;
  //!Returns a QString containing the requested value for the the table in the statistics dialog
  /*!
   * Example :
   *  ____________________________
   * |             |Name   |Cube |
   * |             |_______|_____|
   * |General Info |\#Edges|12   |
   * |_____________|_______|_____|
   * compute stats(0) should return "Cube" and compute_stats(1) should return QString::number(12);
   * The numbers must be coherent with the order of declaration of the titles in the header.
   *
   *
   */
  virtual QString compute_stats(int i);

  //!Contains the number of group and subgroups containing this item.
  int has_group;

public Q_SLOTS:

  //! Notifies the program that the internal data or the properties of
  //! an item has changed, and that it must be computed again.It is
  //! important to call this function whenever the internal data is changed,
  //! or the displayed item will not be updated.
  //!Must be overloaded.
  virtual void invalidateOpenGLBuffers();
  //!Setter for the color of the item. Calls invalidateOpenGLBuffers() so the new color is applied.
  virtual void setColor(QColor c) { color_ = c; invalidateOpenGLBuffers(); }
  //!Setter for the RGB color of the item. Calls setColor(QColor).
  //!@see setColor(QColor c)
  void setRbgColor(int r, int g, int b) { setColor(QColor(r, g, b)); }
  //!Setter for the name of the item.
  virtual void setName(QString n) { name_ = n; }
  //!Setter for the visibility of the item.
  virtual void setVisible(bool b) { visible_ = b; }
  //!Setter for the rendering mode of the item.
  //!@see RenderingMode
  virtual void setRenderingMode(RenderingMode m) { 
    if (supportsRenderingMode(m))
      rendering_mode = m; 
    Q_EMIT redraw();
  }
  //!Set the RenderingMode to Points.
  void setPointsMode() {
    setRenderingMode(Points);
  }
  //!Set the RenderingMode to Wireframe.
  void setWireframeMode() {
    setRenderingMode(Wireframe);
  }

  //!Set the RenderingMode to Flat.
  void setFlatMode() {
    setRenderingMode(Flat);
  }
  //!Set the RenderingMode to FlatPlusEdges.
  void setFlatPlusEdgesMode() {
    setRenderingMode(FlatPlusEdges);
  }
  //!Set the RenderingMode to Gouraud.
  void setGouraudMode() {
    setRenderingMode(Gouraud);
  }
  //!Set the RenderingMode to PointsPlusNormals.
  void setPointsPlusNormalsMode(){
    setRenderingMode(PointsPlusNormals);
  }
  //!Set the RenderingMode to Splatting.
  void setSplattingMode(){
    setRenderingMode(Splatting);
  }

  //! If b is true, the item will use buffers to render the color.
  //! If b is false, it will use a uniform value. For example, when
  //! using the mesh segmentation plugin, the item must be multicolor.
  void setItemIsMulticolor(bool b){
    is_monochrome = !b;
  }
  
  //!Emits an aboutToBeDestroyed() signal.
  virtual void itemAboutToBeDestroyed(Scene_item*);

  //!Selects a point through raycasting.
  virtual void select(double orig_x,
                      double orig_y,
                      double orig_z,
                      double dir_x,
                      double dir_y,
                      double dir_z);

Q_SIGNALS:
  void itemChanged();
  void aboutToBeDestroyed();
  void redraw();

protected:
  //!Holds the BBox of the item
  mutable Bbox _bbox;
  mutable bool is_bbox_computed;
  virtual void compute_bbox()const{}
  // The four basic properties
  //!The name of the item.
  QString name_;
  //!The color of the item.
  QColor color_;
  //!The visibility of the item.
  bool visible_;
  //!Specifies if the item is currently selected.
  bool is_selected;
  //! Specifies if the item is monochrome and uses uniform attribute for its color
  //! or is multicolor and uses buffers.
  bool is_monochrome;
  //! Holds the number of vertices that are not linked to the polyhedron from the OFF
  //! file.
  std::size_t nb_isolated_vertices;
  /*! Decides if the draw function must call initialize_buffers() or not. It is set
   * to true in the end of initialize_buffers() and to false in invalidateOpenGLBuffers(). The need of
   * this boolean comes from the need of a context from the OpenGLFunctions used in
   * initialize_buffers().
   * @see initialize_buffers()
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
   * \brief
   *  Contains the current RenderingMode.
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

  //! Used pass data to the shader.
  int vertexLoc;
  //! Used pass data to the shader.
  int normalLoc;
  //! Used pass data to the shader.
  int colorLoc;
  /*! Fills the VBOs with data. Must be called after each call to #compute_elements().
   * @see compute_elements()
   */
  void initialize_buffers(){}

  /*! Collects all the data for the shaders. Must be called in #invalidateOpenGLBuffers().
   * @see invalidateOpenGLBuffers().
   */
  void compute_elements(){}
  /*! Passes all the uniform data to the shaders.
   * According to program_name, this data may change.
   */
  void attrib_buffers(CGAL::Three::Viewer_interface*, int program_name) const;

  /*! Compatibility function. Calls `viewer->getShaderProgram()`. */
  virtual QOpenGLShaderProgram* getShaderProgram(int name , CGAL::Three::Viewer_interface *viewer = 0) const;
}; // end class Scene_item
}
}

#include <QMetaType>
Q_DECLARE_METATYPE(CGAL::Three::Scene_item*)

#endif // SCENE_ITEM_H
