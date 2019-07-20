#ifndef VIEWER_H
#define VIEWER_H
#include <QMap>
#include <CGAL/Qt/qglviewer.h>
#include <QOpenGLFunctions_2_1>

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

private:
  Scene* m_pScene;
}; // end class Viewer

#endif // VIEWER_H
