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
#include <QTableWidget>

class Viewer : public QGLViewer 
{
  Q_OBJECT

  CGAL::Timer timer;
  Scene* scene;
  bool wireframe;
  bool flatShading;
  bool edges;
  bool vertices;
  unsigned int modeFilledFacet;
  int selectedVolumeIndex;

  typedef LCC::Dart_handle Dart_handle;
  typedef LCC::Dart_const_handle Dart_const_handle;

  std::vector<std::pair<int,Dart_handle> >* pVolumeDartIndex;
  std::vector<char>* pVolumeProperties;


public:
  Viewer(QWidget* parent)
    : QGLViewer(parent), wireframe(false), flatShading(true),
      edges(true), vertices(true), modeFilledFacet(0), selectedVolumeIndex(-1)
  {}

  void setScene(Scene* scene_)
  {
    scene = scene_;
  }

  void setVectorPointers(std::vector<std::pair<int,Dart_handle> >* v1,
                         std::vector<char>* v2)
  {
    pVolumeDartIndex = v1;
    pVolumeProperties = v2;
  }

  void setSelectedVolumeIndex(int index)
  {
    selectedVolumeIndex = index;
  }

public:
  void draw();

  virtual void init();

  void keyPressEvent(QKeyEvent *e);
  
  virtual QString helpString() const;

public slots :
  
  void sceneChanged();

protected:
  void drawFacet(Dart_const_handle ADart);
  void drawEdges(Dart_const_handle ADart);
  void draw_one_vol(Dart_const_handle ADart, bool filled);
  CGAL::Bbox_3 bbox();
};

#endif
