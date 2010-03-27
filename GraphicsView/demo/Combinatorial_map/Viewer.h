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
  unsigned int modeFilledFace;
  int markVolume;
  Map::All_darts_iterator iteratorAllDarts;
  
  typedef Map::Dart_handle Dart_handle;

public:
  Viewer(QWidget* parent)
    : QGLViewer(parent), wireframe(false), flatShading(true),
    edges(true), vertices(true), modeFilledFace(0)
    {}

  void setScene(Scene* scene_)
  {
    scene = scene_;
    markVolume=scene->map.get_new_mark();
    iteratorAllDarts=scene->map.darts_begin();
  }

  void clear();

public:
  void draw();

  virtual void init();
  // void  gl_draw_surface();

  void keyPressEvent(QKeyEvent *e);
  
  virtual QString helpString() const;

public slots :
  
  void sceneChanged();

 protected:
  void drawFace(Dart_handle ADart, int AMark);
  void drawEdges(Dart_handle ADart);

  void draw_one_vol_filled_faces(Dart_handle ADart,
				 int amarkvol, int amarkface);
  
  void draw_current_vol_filled_faces(Dart_handle ADart);
  void draw_current_vol_and_neighboors_filled_faces(Dart_handle ADart);  
};

#endif
