#ifndef VIEWER_H
#define VIEWER_H

#include "Scene.h"
#include <QGLViewer/qglviewer.h>


class Viewer : public QGLViewer {

  typedef qglviewer::Vec Vec;

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

  void init();
  void clear();

public:
  void draw();

  void gl_draw_surface();


public slots :

  void sceneChanged();
  void render_video();
  
signals:
  void valueChanged(int i);

private:
  Vec next_around_circle(const float& phi, const Vec& pos, const Vec& ori);
};

#endif
