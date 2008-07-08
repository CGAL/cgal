#ifndef VIEWER_H
#define VIEWER_H

#include <QGLViewer/qglviewer.h>

// forward declarations
class QWidget;
class Scene;

class Viewer : public QGLViewer {

  Q_OBJECT

public:
  Viewer(QWidget * parent);
  void draw();
  void setScene(Scene* scene);

private:
  Scene* scene;
}; // end class Viewer

#endif // VIEWER_H

