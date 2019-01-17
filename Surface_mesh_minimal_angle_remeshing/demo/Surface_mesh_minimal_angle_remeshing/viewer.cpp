#include "viewer.h"
#include "scene.h"
#include <QMouseEvent>
#include <QGLFunctions>
#include <CGAL/Qt/CreateOpenGLContext.h>

Viewer::Viewer(QWidget *parent)
  : CGAL::QGLViewer(parent),
  m_pScene(NULL),
  m_custom_mouse(false) {
}

void Viewer::setScene(Scene *pScene) {
  this->m_pScene = pScene;
}

void Viewer::draw() {
  CGAL::QGLViewer::draw();
  if (m_pScene != NULL) {
    m_pScene->draw(this);
  }
}

void Viewer::initializeGL() {
  CGAL::QGLViewer::initializeGL();
  setBackgroundColor(::Qt::white);
  //m_pScene->initGL(this);
}

/*void Viewer::mousePressEvent(QMouseEvent *e) {
  if (e->modifiers() == Qt::ControlModifier) {
    m_custom_mouse = true;
  }
  CGAL::QGLViewer::mousePressEvent(e);
}

void Viewer::mouseReleaseEvent(QMouseEvent *e) {
  if (m_custom_mouse) {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QApplication::restoreOverrideCursor();
    m_custom_mouse = false;
  }
  CGAL::QGLViewer::mouseReleaseEvent(e);
}*/
