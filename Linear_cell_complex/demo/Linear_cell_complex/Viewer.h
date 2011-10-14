// Copyright (c) 2010 CNRS, LIRIS, http://liris.cnrs.fr/, All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://gdamiand@scm.gforge.inria.fr/svn/cgal/branches/features/Linear_cell_complex-gdamiand/Linear_cell_complex/demo/Linear_cell_complex/Viewer.h $
// $Id: Viewer.h 58880 2010-09-24 19:41:06Z gdamiand $
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
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
  unsigned int modeFilledFacet;
  int markVolume;
  Map::Dart_range::iterator iteratorAllDarts;
  
  typedef Map::Dart_handle Dart_handle;

public:
 Viewer(QWidget* parent)
    : QGLViewer(parent), wireframe(false), flatShading(true),
    edges(true), vertices(true), modeFilledFacet(0)
    {}

  void setScene(Scene* scene_)
  {
    scene = scene_;
    markVolume=scene->map->get_new_mark();
    iteratorAllDarts=scene->map->darts().begin();
  }

  Map::Dart_range::iterator getCurrentDart() const
    { return iteratorAllDarts; }

  //  void clear();

public:
  void draw();

  virtual void init();
  // void  gl_draw_surface();

  void keyPressEvent(QKeyEvent *e);
  
  virtual QString helpString() const;

public slots :
  
  void sceneChanged();

 protected:
  void drawFacet(Dart_handle ADart, int AMark);
  void drawEdges(Dart_handle ADart);

  void draw_one_vol_filled_facets(Dart_handle ADart,
				 int amarkvol, int amarkfacet);
  
  void draw_current_vol_filled_facets(Dart_handle ADart);
  void draw_current_vol_and_neighboors_filled_facets(Dart_handle ADart);  
};

#endif
