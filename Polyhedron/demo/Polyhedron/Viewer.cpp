#include "Viewer.h"
#include "Scene.h"

Viewer::Viewer(QWidget* parent, bool antialiasing)
  : QGLViewer(parent),
    scene(0),
    antialiasing(antialiasing)
{
  setBackgroundColor(::Qt::white);
}

void Viewer::setScene(Scene* scene)
{
  this->scene = scene;
}

void Viewer::setAntiAliasing(bool b)
{
  antialiasing = b;
  updateGL();
}

void Viewer::draw()
{
  QGLViewer::draw();
  if(scene == 0)
    return;

  ::glLineWidth(1.0f);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.0f,1.0f);
  ::glClearColor(1.0f,1.0f,1.0f,0.0f);
  ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  if(antiAliasing())
  {
    ::glEnable(GL_BLEND);
    ::glEnable(GL_LINE_SMOOTH);
    ::glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }
  else
  {
    ::glDisable(GL_BLEND);
    ::glDisable(GL_LINE_SMOOTH);
    ::glBlendFunc(GL_ONE, GL_ZERO);
    ::glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
  }
  scene->draw();
}
