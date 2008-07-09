#include "Viewer.h"
#include "Scene.h"

Viewer::Viewer(QWidget* parent, bool antialiasing)
  : QGLViewer(parent),
    scene(0),
    antialiasing(antialiasing)
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

  if(antiAliasing())
  {
    ::glEnable(GL_LINE_SMOOTH);
    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    ::glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
  }
  else
  {
    ::glDisable(GL_LINE_SMOOTH);
    ::glDisable (GL_BLEND);
    ::glBlendFunc(GL_ONE, GL_ZERO);
    ::glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
  }
  scene->draw();
}
