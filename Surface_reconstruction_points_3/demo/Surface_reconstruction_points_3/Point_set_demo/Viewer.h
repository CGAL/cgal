#ifndef VIEWER_H
#define VIEWER_H

#include <QGLViewer/qglviewer.h>

// forward declarations
class QWidget;
class Scene_draw_interface;

class Viewer : public QGLViewer {

  Q_OBJECT

public:
  Viewer(QWidget * parent, bool antialiasing = false);

  // overload several QGLViewer virtual functions
  void draw();
  void initializeGL();
  void drawWithNames();
  void postSelection(const QPoint&);

  void setScene(Scene_draw_interface* scene);
  bool antiAliasing() const { return antialiasing; }

signals:
  void selected(int);

public slots:
  void setAntiAliasing(bool b);
  void setTwoSides(bool b);

private:
  void draw_aux(bool with_names);

  Scene_draw_interface* scene;
  bool antialiasing;
  bool twosides;
  bool m_isInitialized;

}; // end class Viewer

#endif // VIEWER_H
