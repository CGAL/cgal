#ifndef VIEWER_H
#define VIEWER_H

#include "typedefs.h"
#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>

class Viewer : public QGLViewer {
  Q_OBJECT

  CGAL::Timer timer;
  Scene* scene;
  bool wireframe;
  bool flatShading;
  bool edges;
  bool vertices;

  typedef Map::Dart_handle Dart_handle;

public:
  Viewer(QWidget* parent)
    : QGLViewer(parent), wireframe(false), flatShading(true),
    edges(true), vertices(true)
  {}

  void setScene(Scene* scene_)
  { scene = scene_; }

  void clear();

public:
  void draw();
  void drawFace(Dart_handle ADart, int AMark);
  void drawEdges(Dart_handle ADart);

  virtual void init();
  // void  gl_draw_surface();

  void keyPressEvent(QKeyEvent *e);
  
  virtual QString helpString() const;

public slots :
  
  void sceneChanged();   
};

#endif
