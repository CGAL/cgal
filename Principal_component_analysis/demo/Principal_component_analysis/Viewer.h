#ifndef VIEWER_H
#define VIEWER_H
#include <QMap>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLFunctions_2_1>

// forward declarations
class QWidget;
class Scene;

class Viewer : public QGLViewer,
public QOpenGLFunctions_2_1{

  Q_OBJECT

public:
  Viewer(QWidget * parent);

  // overload several QGLViewer virtual functions
  void draw();
  void initializeGL();
  void setScene(Scene* pScene);

private:
  Scene* m_pScene;
}; // end class Viewer

#endif // VIEWER_H
