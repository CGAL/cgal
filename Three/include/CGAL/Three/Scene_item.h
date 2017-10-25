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
#include <QString>
#include <QPixmap>
#include <vector>
#include <CGAL/Three/Scene_interface.h>  //gives access to RenderingMode
#include <CGAL/Bbox_3.h>

struct D;
class QOpenGLFramebufferObject;
namespace CGAL {
namespace Three {
  class Viewer_interface;
  struct Triangle_container;
  struct Edge_container;
}
}
namespace qglviewer {
  class ManipulatedFrame;
}

class QMenu;
class QKeyEvent;
class QSlider;
namespace CGAL {
namespace Three {

class Scene_group_item;
class Viewer_interface;
//! This class represents an object in the OpenGL scene.
//! It contains all the functions called by the Scene.
class SCENE_ITEM_EXPORT Scene_item : public QObject{
  Q_OBJECT
  Q_PROPERTY(QColor color READ color WRITE setColor)
  Q_PROPERTY(QString name READ name WRITE setName)
  Q_PROPERTY(bool visible READ visible WRITE setVisible)
  Q_ENUMS(RenderingMode)
  Q_PROPERTY(RenderingMode renderingMode READ renderingMode WRITE setRenderingMode)
public:
  typedef CGAL::Bbox_3 Bbox;
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;
  //!
  //! \brief The Gl_data_name enum is used as a flag to specify what should be
  //! re-computed during `computeElements()`. The flag correspondig to this enum is
  //! `Gl_data_names`, and multiple flags can be combined whith the operator |.
  //! For instance, you can use GEOMETRY|COLORS as a single value.
  //!
  enum Gl_data_name{
    GEOMETRY = 0x1,                     //!< Invalidates the vertices, edges and faces.
    COLORS   = 0x2,                     //!< Invalidates the color of each vertex
    NORMALS  = 0x4,                     //!< Invalidate the normal of each vertex.
    ALL      = GEOMETRY|COLORS|NORMALS  //!< Invalidate everything
  };

#ifdef DOXYGEN_RUNNING
  //! \brief Flag interface for Scene_item::Gl_data_name.
  enum Gl_data_names{};
#endif
  Q_DECLARE_FLAGS(Gl_data_names, Gl_data_name)

  //! \brief The default color of a scene_item.
  //!
  //! This color is the one that will be displayed if none is specified after its creation.
  static const QColor defaultColor();
  //!
  //! \brief selectionColor is the color used for selected item.
  //!
  QColor selectionColor();
  Scene_item();

  virtual ~Scene_item();
  //! \brief Duplicates the item.
  //!
  //! Creates a new item as a copy of this one.
  //! Must be overriden;
  virtual Scene_item* clone() const = 0;

  //! \brief Indicates if `m` is supported
  //!
  //! If it is, it will be displayed in the context menu of the item.
  //! Must be overriden;
  virtual bool supportsRenderingMode(RenderingMode m) const = 0;
  /*! \brief The drawing function for faces.
   *
   * Draws the faces of the item in the viewer.
   * \param viewer the active `Viewer_interface`
   * \param pass the current pass in the Depth Peeling (transparency) algorithm.
   * -1 means that no depth peeling is applied.
   * \param writing_depth means that the color of the faces will be drawn in a grayscale
   * according to the depth of the fragment in the shader. It is used by the transparency.
   * \param fbo contains the texture used by the Depth Peeling algorithm.
   * Should be NULL if pass <= 0;
   */  
#ifdef DOXYGEN_RUNNING
    virtual void draw(CGAL::Three::Viewer_interface* viewer,
                      int pass, bool writing_depth, QOpenGLFramebufferObject* fbo)
#else
  virtual void draw(CGAL::Three::Viewer_interface* ,
                    int , bool , QOpenGLFramebufferObject* )
#endif
  {}
  /*! \brief The drawing function for the edges.
   *
   * Draws the edges and lines of the item in the viewer.
   */
  virtual void drawEdges(CGAL::Three::Viewer_interface* ) {}
  /*! \brief The drawing function for the points.
   *
   * Draws the points of the item in the viewer.
   */
  virtual void drawPoints(CGAL::Three::Viewer_interface* ) {}

  //! Called by the scene. If b is true, then this item is currently selected.
  virtual void selection_changed(bool b);

  // Functions for displaying meta-data of the item
  //!\brief Contains meta-data about the item.
  //!
  //! Must be overriden;
  //! @returns a QString containing meta-data about the item.
  virtual QString toolTip() const = 0;
  //! \brief Contains graphical meta-data about the item.
  //! @returns a QPixmap containing graphical meta-data about the item.
  virtual QPixmap graphicalToolTip() const { return QPixmap(); }

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
  virtual Bbox bbox() const = 0 ;


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
  virtual QColor color() const;
  //!Getter for the item's name.
  //! @returns the current name of the item.
  virtual QString name() const;
  //! If the item is not visible, it is not drawn.
  //! @returns the current visibility of the item.
  virtual bool visible() const;
  //!Getter for the item's rendering mode.
  //! @returns the current rendering mode of the item.
  //!@see RenderingMode
  virtual RenderingMode renderingMode() const;
  //!The renderingMode's name.
  //! @returns the current rendering mode of the item as a human readable string.
  virtual QString renderingModeName() const;

  //! \brief The context menu of an item.
  //!
  //! Contains the list of the supported rendering modes,
  //! the Operations menu, actions to save or clone the item if it is supported
  //! and any contextual action for the item.
  virtual QMenu* contextMenu();
  //!
  //! \brief isWriting is a property associated with threads.
  //! \return `true` if the item is being modified.
  //! \see `writing()`
  //!
  bool isWriting()const;
  //!
  //! \brief writing sets the state of the item to being modified.
  //!
  //! Thread safe.
  //! Some features will be disabled.
  //!\see `doneWriting()`
  //! \see `isWriting()`
  void writing();
  //!
  //! \brief doneWriting sets the state of the item to done being modified.
  //!
  //! Thread safe.
  //!
  void doneWriting();
  //!
  //! \brief isReading is a property associated with threads.
  //!
  //! Thread safe.
  //! \return `true` if something is parsing the item data.
  //! \see `reading()`
  //!
  int isReading()const;
  //!
  //! \brief reading sets the state of the item to being parsed.
  //!
  //! Thread safe.
  //! Some features will be disabled.
  //!
  void reading();
  //!
  //! \brief doneReading sets the state of the item to done being parsed.
  //!
  //! Thread safe.
  //!
  void doneReading();

  //!
  //! \brief setId informs the item of its current index in the scene entries.
  //!
  void setId(int id);

  //!
  //! \brief getId returns the current index of this item in the scene entries.
  //!
  int getId()const;

  //!Handles key press events received from the event filters installed for this item.
  virtual bool keyPressEvent(QKeyEvent*){return false;}

  //!The group containing the item.
  //! \returns the parent group if the item is in a group
  //! \returns 0 if the item is not in a group.
  Scene_group_item* parentGroup() const;

  //!Contains the header for the table in the statistics dialog
  /*!
     A header data is composed of 2 columns : the Categories and the titles.
     A category is the name given to an association of titles.
     A title is the name of a line.
     \verbatim
     For example,
     Category :    | Titles| Values
     2 lines       |       |
      ____________________________
     |             |Name   |Cube |
     |             |_______|_____|
     |General Info | #Edges|12   |
     |_____________|_______|_____|

      would be stored as follows :
     categories = std::pair<QString,int>(QString("General Info"),2)
     titles.append("Name");
     titles.append("#Edges");\endverbatim
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
  int hasGroup() const;
  //!Sets the number of group and subgroups containing this item.
  void setHasGroup(int );
  //!
  //! \brief newViewer adds Vaos for `viewer`.
  //!
  //! Must be overriden;
  //!
  virtual void newViewer(CGAL::Three::Viewer_interface* viewer) = 0;
  //!
  //! \brief removeViewer removes the Vaos fo `viewer`.
  //!
  //! Must be overriden;
  //!
  virtual void removeViewer(CGAL::Three::Viewer_interface* viewer) = 0;

  //! Returns the selection status of this item.
  bool isSelected() const;

  //!Set the selection status of this item
  void setSelected(bool b);

  /*! Collects all the data for the shaders. Must be called in #invalidate().
   * @see invalidate().
   */
  virtual void computeElements(Gl_data_names){}

public Q_SLOTS:

  //! Notifies the application that the internal data corresponding to `name` is not valid anymore,
  //! and asks for a re-computation of this data. It is
  //! important to call this function whenever the internal data is changed,
  //! or the displayed item will not be updated.
  virtual void invalidate(Gl_data_names name);
  //!Setter for the color of the item.
  virtual void setColor(QColor c);
  //!Setter for the RGB color of the item. Calls setColor(QColor).
  //!@see setColor(QColor c)
  void setRbgColor(int r, int g, int b){ setColor(QColor(r, g, b)); }
  //!Sets the name of the item.
  virtual void setName(QString n);
    //!Sets the visibility of the item.
  virtual void setVisible(bool b);
  //!Set the parent group. If `group==0`, then the item has no parent.
  //!This function is called by `Scene::changeGroup` and should not be
  //!called manually.
  virtual void moveToGroup(Scene_group_item* group);
  //!Sets the rendering mode of the item.
  //!@see RenderingMode
  virtual void setRenderingMode(RenderingMode m) ;
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

  //!Returns the alpha value for the item.
  //! Must be called within a valid openGl context.
  virtual float alpha() const;

  //! Sets the value of the aplha Slider for this item.
  //!
  //! Must be overriden;
  //! \param alpha must be between 0 and 255
  virtual void setAlpha(int alpha) = 0;

  //!Selects a point through raycasting.
  virtual void select(double orig_x,
                      double orig_y,
                      double orig_z,
                      double dir_x,
                      double dir_y,
                      double dir_z);



Q_SIGNALS:
  //! Is emitted to notify a change in the item's geometric data.
  //! Will replace the information of this item and re-draw the scene.
  void itemChanged();
  //! Is emitted when the item is shown to notify a change in the item's visibility.
  //! Typically used to update the scene's bbox;
  void itemVisibilityChanged();
  //! Is emitted to notify that the item is about to be deleted.
  void aboutToBeDestroyed();
  //! Is emitted to require a new display.
  //! Use it to avoid having to move the camera to see the changes of this item.
  void redraw();
  //!
  //! Is emitted when the data computation of the item is finished and the thread
  //! that were performing it is done.
  //!
  void dataProcessed();

protected:
private:
  friend struct D;
  mutable D* d;
}; // end class Scene_item
Q_DECLARE_OPERATORS_FOR_FLAGS(CGAL::Three::Scene_item::Gl_data_names)
}
}

#include <QMetaType>
Q_DECLARE_METATYPE(CGAL::Three::Scene_item*)
#endif // SCENE_ITEM_H

