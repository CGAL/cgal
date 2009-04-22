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
  draw_aux(false);
}

void Viewer::initializeGL()
{
  QGLViewer::initializeGL();
  scene->initializeGL();
}

void Viewer::draw_aux(bool with_names)
{
  QGLViewer::draw();
  if(scene == 0)
    return;

  ::glLineWidth(1.0f);
  ::glPointSize(10.f);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.0f,1.0f);
  ::glClearColor(1.0f,1.0f,1.0f,0.0f);
  ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  if(antiAliasing())
  {
    ::glEnable(GL_BLEND);
    ::glEnable(GL_LINE_SMOOTH);
    ::glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }
  else
  {
    ::glDisable(GL_BLEND);
    ::glDisable(GL_LINE_SMOOTH);
    ::glDisable(GL_POLYGON_SMOOTH_HINT);
    ::glBlendFunc(GL_ONE, GL_ZERO);
    ::glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
  }
  scene->draw(with_names);
}

void Viewer::drawWithNames()
{
  draw_aux(true);
}

void Viewer::postSelection(const QPoint&)
{
  emit selected(this->selectedName());
}
