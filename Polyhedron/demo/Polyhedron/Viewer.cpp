#include "Viewer.h"
#include "Scene.h"

Viewer::Viewer(QWidget* parent)
  : QGLViewer(parent), scene(0)
{
  setBackgroundColor(Qt::white);
}

void 
Viewer::setScene(Scene* scene)
{
  this->scene = scene;
}

void
Viewer::draw()
{
  QGLViewer::draw();
  if(scene == 0)
    return;

  ::glClearColor(1., 1., 1., 0.);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.0f,1.0f);
  ::glLineWidth(1.0f);
  ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  ::glEnable(GL_LINE_SMOOTH);
  scene->draw();
}
