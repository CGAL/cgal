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
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Kumar Snehasish <kumar.snehasish@gmail.com>
//
#ifndef VIEWER_H
#define VIEWER_H

#include "typedefs.h"
#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>

class Viewer : public QGLViewer
{
  Q_OBJECT

  CGAL::Timer timer;
  Scene* scene;
  bool wireframe;
  bool flatShading;
  bool edges;
  bool vertices;
  CGAL::Bbox_3 bb;

  GLuint m_dlFaces;
  GLuint m_dlFacesFlat;
  GLuint m_dlEdges;
  GLuint m_dlVertices;
  bool m_displayListCreated;

  typedef LCC::Dart_handle Dart_handle;
  typedef LCC::Dart_const_handle Dart_const_handle;


public:
  Viewer(QWidget* parent)
    : QGLViewer(parent), wireframe(false), flatShading(true),
      edges(true), vertices(true), m_displayListCreated(false)
  {
    QGLFormat newFormat = this->format();
    newFormat.setSampleBuffers(true);
    newFormat.setSamples(16);
    this->setFormat(newFormat);
  }

  void setScene(Scene* scene_)
  {
    scene = scene_;
  }

public:
  void draw();

  virtual void init();

  void keyPressEvent(QKeyEvent *e);

  virtual QString helpString() const;

public slots :

  void sceneChanged();

protected:
  void initDraw();
  void drawAllFaces(bool flat);
  void drawAllEdges();
  void drawAllVertices();
};

#endif
