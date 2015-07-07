#ifndef VIEWER_H
#define VIEWER_H

#include <QGLViewer/qglviewer.h>

// forward declarations
class QWidget;
class Scene;

class Viewer : public QGLViewer{

  Q_OBJECT

public:
  Viewer(QWidget * parent);

  // overload several QGLViewer virtual functions
  void draw();
  void initializeGL();
  void setScene(Scene* pScene);
  QOpenGLContext *oglContext()const {return oglContext_;}

protected:
  virtual void mousePressEvent(QMouseEvent* e);
  virtual void mouseReleaseEvent(QMouseEvent* e);
  
private:
  static QGLContext* createContext()
  {
      QOpenGLContext *context = new QOpenGLContext();
      QSurfaceFormat format;
      format.setVersion(3,3);
      format.setProfile(QSurfaceFormat::CompatibilityProfile);
      context->setFormat(format);
      return QGLContext::fromOpenGLContext(context);
  }

  Scene* m_pScene;
  bool m_custom_mouse;
  QOpenGLContext *oglContext_;
}; // end class Viewer

#endif // VIEWER_H
