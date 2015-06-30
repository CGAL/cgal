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

/* TODO: code to draw non convex faces (since glu tesselator is now deprecated)

  #include <GL/glu.h>
  #include <GL/glut.h>
  #ifndef CALLBACK
  #define CALLBACK
  #endif

static void CALLBACK monCombine(GLdouble c[3], void * d [4], GLfloat w [4], void **out)
{
   GLdouble *nv = (GLdouble *) malloc(sizeof(GLdouble)*3);
   nv[0] = c[0];
   nv[1] = c[1];
   nv[2] = c[2];
   *out = nv;
}
*/

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
  bool m_previous_scene_empty;

  //  GLUtesselator* FTess;

  typedef LCC::Dart_handle Dart_handle;
  typedef LCC::Dart_const_handle Dart_const_handle;

public:
  Viewer(QWidget* parent)
    : QGLViewer(parent), wireframe(false), flatShading(true),
      edges(true), vertices(true), m_displayListCreated(false),
      m_previous_scene_empty(true)
  {
    QGLFormat newFormat = this->format();
    newFormat.setSampleBuffers(true);
    newFormat.setSamples(16);
    this->setFormat(newFormat);

    /*    FTess = gluNewTess();
    gluTessCallback(FTess, GLU_TESS_BEGIN,   (void (CALLBACK*)()) &glBegin);
    gluTessCallback(FTess, GLU_TESS_END,     (void (CALLBACK*)()) &glEnd);
    gluTessCallback(FTess, GLU_TESS_VERTEX,  (void (CALLBACK*)()) &glVertex3dv);
    gluTessCallback(FTess, GLU_TESS_COMBINE, (void (CALLBACK*)()) &monCombine); */
  }

  ~Viewer()
  {
    // gluDeleteTess(FTess);
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

public Q_SLOTS:

  void sceneChanged();

protected:
  void initDraw();
  void drawAllFaces(bool flat);
  void drawAllEdges();
  void drawAllVertices();
  void drawOneFaceWireframe(Dart_handle);
  void drawOneFilledFace(Dart_handle, bool flat);
};

#endif
