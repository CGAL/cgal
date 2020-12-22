// Copyright (c) 2009,2010,2012,2015  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

#ifndef SCENE_DRAW_INTERFACE_H
#define SCENE_DRAW_INTERFACE_H

#include <CGAL/license/Three.h>

#include <QPoint>
#include <QVector3D>
class QKeyEvent;
class QPoint;
namespace CGAL
{
namespace Three {
  class Viewer_interface;
  class Scene_item;


//! Base class to interact with the scene from the viewer.
class Scene_draw_interface {
public:
  virtual ~Scene_draw_interface(){}

  /*! Is called by Viewer::initializeGL(). Allows all the initialization
   * of OpenGL code that needs a context.
   */
  virtual void initializeGL(CGAL::Three::Viewer_interface*) = 0;

  //! \brief Draws the items.
  //! It is called by Viewer::draw().
  virtual void draw(CGAL::Three::Viewer_interface*) = 0;
  //!\brief draws the scene in a hidden frame to perform picking.
  //! Is called by Viewer::drawWithNames().
  virtual void drawWithNames(CGAL::Three::Viewer_interface*) = 0;
  //!Pick the point `e` on the screen.
  virtual void setPickedPixel(const QPoint &e) = 0;
  //! \brief Manages the key events.
  //! Override this function to perform actions when keys are pressed.
  //! @returns true if the keyEvent executed well.
  //!
  virtual bool keyPressEvent(QKeyEvent* e) = 0;
  //!\brief print theTextItems.
  virtual void printPrimitiveId(QPoint point, CGAL::Three::Viewer_interface*) = 0;
  //!\brief update theTextItems.
  virtual void updatePrimitiveIds(CGAL::Three::Scene_item*) = 0;

  /*!
   * \brief checks if the text at position (x,y,z) is visible or not.
   * \param x the X coordinate of theTextItem's position.
   * \param y the Y coordinate of theTextItem's position.
   * \param z the Z coordinate of theTextItem's position.
   * \param viewer the viewer used to display the Scene.
   * \return true if the TextItem is visible. */
  virtual bool  testDisplayId(double x, double y, double z, CGAL::Three::Viewer_interface* viewer) = 0;

  ///\brief displays all the vertices ids if there are less than max_textItems.
  virtual void printVertexIds() = 0;
  ///\brief displays all the edges ids if there are less than max_textItems.
  virtual void printEdgeIds() = 0;
  ///\brief displays all the faces ids if there are less than max_textItems.
  virtual void printFaceIds() = 0;
  ///\brief displays all the primitive ids if there are less than max_textItems.
  virtual void printAllIds() = 0;

  //!\brief moves the camera orthogonally to the picked face.
  //!
  //! \param point the picked point
  //! \param viewer the active viewer
  virtual void zoomToPosition(QPoint point,
                        CGAL::Three::Viewer_interface* viewer) = 0;
};
}
}
#endif // SCENE_DRAW_INTERFACE_H;
