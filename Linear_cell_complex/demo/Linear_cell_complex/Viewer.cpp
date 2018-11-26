// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Kumar Snehasish <kumar.snehasish@gmail.com>
//
#include "Viewer.h"

#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Qt/CreateOpenGLContext.h>
#include <CGAL/Qt/viewer_actions.h>
#include <CGAL/Qt/vec.h>
#include <QDebug>

Viewer::Viewer(QWidget* parent) :
    Base(parent, NULL, ""),
    wireframe(false),
    flatShading(true),
    edges(true),
    vertices(true),
    inverse_normal(false),
    size_points(7.),
    size_edges(3.1),
    ambient(0.6f, 0.5f, 0.5f, 0.5f),
    m_previous_scene_empty(true)
{}

Viewer::~Viewer()
{}

void Viewer::sceneChanged()
{
  Base::compute_elements();
  this->camera()->setSceneBoundingBox(CGAL::qglviewer::Vec(m_bounding_box.xmin(),
                                                     m_bounding_box.ymin(),
                                                     m_bounding_box.zmin()),
                                      CGAL::qglviewer::Vec(m_bounding_box.xmax(),
                                                     m_bounding_box.ymax(),
                                                     m_bounding_box.zmax()));
  if (m_previous_scene_empty)
    this->showEntireScene();
  else
    this->update();

  m_previous_scene_empty = scene->lcc->is_empty(); // for the next call to sceneChanged
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
  const Qt::KeyboardModifiers modifiers = e->modifiers();

  if ((e->key()==Qt::Key_W) && (modifiers==Qt::NoButton))
  {
    wireframe = !wireframe;
    if (wireframe)
    {
      displayMessage("Wireframe.");
    }
    else
    {
      displayMessage("Filled faces.");
    }
    update();
  }
  else if ((e->key()==Qt::Key_F) && (modifiers==Qt::NoButton))
  {
    flatShading = !flatShading;
    if (flatShading)
      displayMessage("Flat shading.");
    else
      displayMessage("Gouraud shading.");

    update();

  }
  else if ((e->key()==Qt::Key_E) && (modifiers==Qt::NoButton))
  {
    edges = !edges;
    displayMessage(QString("Draw edges=%1.").arg(edges?"true":"false"));

    update();
  }
  else if ((e->key()==Qt::Key_V) && (modifiers==Qt::NoButton))
  {
    vertices = !vertices;
    displayMessage(QString("Draw vertices=%1.").arg(vertices?"true":"false"));
    update();
  }
  else if ((e->key()==Qt::Key_N) && (modifiers==Qt::NoButton))
  {
    inverse_normal = !inverse_normal;
    displayMessage(QString("Inverse normal=%1.").arg(inverse_normal?"true":"false"));
    sceneChanged();
  }
  else if ((e->key()==Qt::Key_Plus) && (modifiers==Qt::KeypadModifier))
  {
    size_edges+=.5;
    displayMessage(QString("Size of edges=%1.").arg(size_edges));
    update();
  }
  else if ((e->key()==Qt::Key_Minus) && (modifiers==Qt::KeypadModifier))
  {
    if (size_edges>.5) size_edges-=.5;
    displayMessage(QString("Size of edges=%1.").arg(size_edges));
    update();
  }
  else if ((e->key()==Qt::Key_Plus) && (modifiers==(Qt::ShiftModifier|Qt::KeypadModifier)))
  {
    size_points+=.5;
    displayMessage(QString("Size of points=%1.").arg(size_points));
    update();
  }
  else if ((e->key()==Qt::Key_Minus) && (modifiers==(Qt::ShiftModifier|Qt::KeypadModifier)))
  {
    if (size_points>.5) size_points-=.5;
    displayMessage(QString("Size of points=%1.").arg(size_points));
    update();
  }
  else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::NoButton))
  {
    ambient.setX(ambient.x()+.1);
    if (ambient.x()>1.) ambient.setX(1.);
    ambient.setY(ambient.x()+.1);
    if (ambient.y()>1.) ambient.setY(1.);
    ambient.setZ(ambient.x()+.1);
    if (ambient.z()>1.) ambient.setZ(1.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
    update();
  }
  else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::NoButton))
  {
    ambient.setX(ambient.x()-.1);
    if (ambient.x()<0.) ambient.setX(0.);
    ambient.setY(ambient.y()-.1);
    if (ambient.y()<0.) ambient.setY(0.);
    ambient.setZ(ambient.z()-.1);
    if (ambient.z()<0.) ambient.setZ(0.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
    update();
  }
  else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::ShiftModifier))
  {
    ambient.setX(ambient.x()+.1);
    if (ambient.x()>1.) ambient.setX(1.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
    update();
  }
  else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::AltModifier))
  {
    ambient.setY(ambient.y()+.1);
    if (ambient.y()>1.) ambient.setY(1.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
    update();
  }
  else if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::ControlModifier))
  {
    ambient.setZ(ambient.z()+.1);
    if (ambient.z()>1.) ambient.setZ(1.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
    update();
  }
  else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::ShiftModifier))
  {
    ambient.setX(ambient.x()-.1);
    if (ambient.x()<0.) ambient.setX(0.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
    update();
  }
  else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::AltModifier))
  {
    ambient.setY(ambient.y()-.1);
    if (ambient.y()<0.) ambient.setY(0.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
    update();
  }
  else if ((e->key()==Qt::Key_PageDown) && (modifiers==Qt::ControlModifier))
  {
    ambient.setZ(ambient.z()-.1);
    if (ambient.z()<0.) ambient.setZ(0.);
    displayMessage(QString("Light color=(%1 %2 %3).").
        arg(ambient.x()).arg(ambient.y()).arg(ambient.z()));
    update();
  }
  else
    CGAL::QGLViewer::keyPressEvent(e);
}

QString Viewer::helpString() const
{
  QString text("<h2>L C C   V i e w e r</h2>");
  text += "Use the mouse to move the camera around the object. ";
  text += "You can respectively revolve around, zoom and translate with "
    "the three mouse buttons. ";
  text += "Left and middle buttons pressed together rotate around the "
    "camera view direction axis<br><br>";
  text += "Pressing <b>Alt</b> and one of the function keys "
    "(<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
  text += "Simply press the function key again to restore it. Several "
    "keyFrames define a ";
  text += "camera path. Paths are saved when you quit the application and "
    "restored at next start.<br><br>";
  text += "Press <b>F</b> to display the frame rate, <b>A</b> for the "
    "world axis, ";
  text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to "
    "save a snapshot. ";
  text += "See the <b>Keyboard</b> tab in this window for a complete "
    "shortcut list.<br><br>";
  text += "Double clicks automates single click actions: A left button "
    "double click aligns the closer axis with the camera (if close enough). ";
  text += "A middle button double click fits the zoom of the camera and "
    "the right button re-centers the scene.<br><br>";
  text += "A left button double click while holding right button pressed "
    "defines the camera <i>Revolve Around Point</i>. ";
  text += "See the <b>Mouse</b> tab and the documentation web pages for "
    "details.<br><br>";
  text += "Press <b>Escape</b> to exit the viewer.";
  return text;
}

