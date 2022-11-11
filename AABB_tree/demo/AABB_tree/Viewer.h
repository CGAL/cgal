#ifndef VIEWER_H
#define VIEWER_H
#include <QMap>
#include <CGAL/Qt/qglviewer.h>


// forward declarations
class QWidget;
class Scene;
class Viewer : public CGAL::QGLViewer{

  Q_OBJECT

public:
  Viewer(QWidget * parent);

  // overload several CGAL::QGLViewer virtual functions
  void draw();
  void initializeGL();
  void setScene(Scene* pScene);

protected:
  virtual void mousePressEvent(QMouseEvent* e);
  virtual void mouseReleaseEvent(QMouseEvent* e);

private:
  Scene* m_pScene;
  bool m_custom_mouse;
}; // end class Viewer

#endif // VIEWER_H
