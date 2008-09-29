#ifndef VIEWER_H
#define VIEWER_H

#include "typedefs.h"
#include <QGLViewer/qglviewer.h>


class Viewer : public QGLViewer {
  Q_OBJECT

  CGAL::Timer timer;
  Scene* scene;

  int nr_of_facets;
public:
  Viewer(QWidget* parent)
    : QGLViewer(parent)
  {}

  void setScene(Scene* scene_)
  {
    scene = scene_;
  }

  void clear();

public:
  void draw();

  void 
  gl_draw_surface();

 public slots :

 void sceneChanged();
   
};

#endif
