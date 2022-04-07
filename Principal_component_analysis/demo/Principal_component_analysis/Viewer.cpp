#include "Viewer.h"
#include "Scene.h"

Viewer::Viewer(QWidget* parent)
  : CGAL::QGLViewer(parent),
    m_pScene(nullptr)
{
}

void Viewer::setScene(Scene* pScene)
{
  this->m_pScene = pScene;
}

void Viewer::draw()
{
  CGAL::QGLViewer::draw();
  if(m_pScene != nullptr)
  {
        glClearColor(1.0f,1.0f,1.0f,1.0f);
        m_pScene->draw(this);
  }
}

void Viewer::initializeGL()
{
  CGAL::QGLViewer::initializeGL();
  //makeCurrent();
  //initializeOpenGLFunctions();
  setBackgroundColor(::Qt::white);
}

