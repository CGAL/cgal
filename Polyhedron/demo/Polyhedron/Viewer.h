#ifndef VIEWER_H
#define VIEWER_H

#include <QGLViewer/qglviewer.h>

// forward declarations
class QWidget;
class Scene;

class Viewer : public QGLViewer {

  Q_OBJECT

public:
  Viewer(QWidget * parent, bool antialiasing = false);
  void draw();
  void setScene(Scene* scene);
  bool antiAliasing() const { return antialiasing; }

public slots:
  void setAntiAliasing(bool b);

private:
  Scene* scene;
  bool antialiasing;
}; // end class Viewer

#endif // VIEWER_H
