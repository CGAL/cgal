#include "Viewer.h"
#include "Scene.h"
#include <QMouseEvent>
#include <QGLFunctions>
#include <CGAL/Qt/CreateOpenGLContext.h>

Viewer::Viewer(QWidget* parent)
  : QGLViewer(CGAL::Qt::createOpenGLContext(),parent),
    m_pScene(NULL),
    m_custom_mouse(false)
{
}

void Viewer::setScene(Scene* pScene)
{
  this->m_pScene = pScene;
}

void Viewer::draw()
{
  glEnable(GL_DEPTH_TEST);
  QGLViewer::draw();
  if(m_pScene != NULL)
  {
      m_pScene->draw(this);
  }

}

void Viewer::initializeGL()
{
  QGLViewer::initializeGL();
  setBackgroundColor(::Qt::white);
  m_pScene->initGL(this);
}

void Viewer::mousePressEvent(QMouseEvent* e)
{
  if ( e->modifiers() == Qt::ControlModifier )
  {
    m_pScene->set_fast_distance(true);
    // Refresh distance function
    m_pScene->cutting_plane(true);
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
    m_pScene->cutting_plane(true);
    QApplication::restoreOverrideCursor();
      
    m_custom_mouse = false;
  }
  
  QGLViewer::mouseReleaseEvent(e);
}

