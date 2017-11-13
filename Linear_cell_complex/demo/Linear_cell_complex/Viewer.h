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

#include <vector>
#include <QMap>
#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>
#include <QOpenGLFunctions_2_1>
#include <QOpenGLVertexArrayObject>
#include <QGLBuffer>
#include <QOpenGLShaderProgram>

#define NB_VBO_BUFFERS 8
#define NB_VAO_BUFFERS 4

class Viewer : public QGLViewer, public QOpenGLFunctions_2_1
{
  Q_OBJECT

  typedef LCC::Dart_handle Dart_handle;
  typedef LCC::Dart_const_handle Dart_const_handle;

public:
  Viewer(QWidget* parent);

  ~Viewer();

  void setScene(Scene* scene_)
  { scene = scene_; }

public:
  void draw();

  virtual void init();

  void keyPressEvent(QKeyEvent *e);

  virtual QString helpString() const;

public Q_SLOTS:

  void sceneChanged();

private:
  void initialize_buffers();
  void attrib_buffers(QGLViewer*);
  void compile_shaders();

  void compute_elements();
  void compute_face(Dart_handle dh, LCC::size_type markface);
  void compute_edge(Dart_handle dh, LCC::size_type markedge);
  void compute_vertex(Dart_handle dh, LCC::size_type markvertex, bool& empty);

private:
  Scene* scene;
  
  bool wireframe;
  bool flatShading;
  bool edges;
  bool vertices;  
  bool inverse_normal;
  
  double size_points;
  double size_edges;

  QVector4D ambient;
  
  bool m_previous_scene_empty;
  bool are_buffers_initialized;

  //Shaders elements
  int vertexLocation[3];
  int normalsLocation;
  int mvpLocation[2];
  int mvLocation;
  int colorLocation;
  int colorsLocation;
  int lightLocation[5];

  std::vector<float> pos_points;
  std::vector<float> pos_lines;
  std::vector<float> pos_facets;
  std::vector<float> smooth_normals;
  std::vector<float> flat_normals;
  std::vector<float> colors;

  QGLBuffer buffers[NB_VBO_BUFFERS];
  QOpenGLVertexArrayObject vao[NB_VAO_BUFFERS];
  QOpenGLShaderProgram rendering_program;
  QOpenGLShaderProgram rendering_program_p_l;

  CGAL::Bbox_3 bb;
};

#endif
