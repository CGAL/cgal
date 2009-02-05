#ifndef VIEWER_H
#define VIEWER_H

#include "typedefs.h"
#include <QGLViewer/qglviewer.h>


class Viewer : public QGLViewer {
  Q_OBJECT

  Scene* scene;
  bool frame_has_been_spun;

public:
  Viewer(QWidget* parent);
  ~Viewer();

  void setScene(Scene* scene_)
  {
    scene = scene_;
  }

  void clear();

public:
  void draw();

  void gl_draw_vertices();

  void gl_draw_surface();

  void gl_draw_constraints();

public slots :
  void frameSpun() {
    frame_has_been_spun = true;
    updateGL();
  }

  void sceneChanged();
   
};

#endif
