#include "Viewer.h"
#include "Scene.h"

Viewer::Viewer(QWidget* parent)
  : QGLViewer(parent),
    m_pScene(NULL)
{
  setBackgroundColor(::Qt::white);
}

void Viewer::setScene(Scene* pScene)
{
  this->m_pScene = pScene;
}

void Viewer::draw()
{
  QGLViewer::draw();
  if(m_pScene != NULL)
  {
	::glClearColor(1.0f,1.0f,1.0f,0.0f);
	m_pScene->draw();
  }
}

void Viewer::initializeGL()
{
  QGLViewer::initializeGL();
}

