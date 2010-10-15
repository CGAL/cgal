#include "Viewer.h"
#include "Scene.h"
#include <QMouseEvent>

Viewer::Viewer(QWidget* parent)
  : QGLViewer(parent),
    m_pScene(NULL),
    m_custom_mouse(false)
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

void Viewer::mousePressEvent(QMouseEvent* e)
{
  if ( e->modifiers() == Qt::ControlModifier )
  {
    m_pScene->set_fast_distance(true);
    // Refresh distance function
    m_pScene->cutting_plane();
    m_custom_mouse = true;
  }
  
  QGLViewer::mousePressEvent(e);
}

void Viewer::mouseReleaseEvent(QMouseEvent* e)
{
  if ( m_custom_mouse )
  {
    m_pScene->set_fast_distance(false);
    
    // Recompute distance function
    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_pScene->cutting_plane();
    QApplication::restoreOverrideCursor();
      
    m_custom_mouse = false;
  }
  
  QGLViewer::mouseReleaseEvent(e);
}


