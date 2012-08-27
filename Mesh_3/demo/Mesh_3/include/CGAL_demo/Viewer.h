#ifndef VIEWER_H
#define VIEWER_H

#include <CGAL_demo/Viewer_config.h>
#include <QGLViewer/qglviewer.h>
#include <QPoint>

// forward declarations
class QWidget;
class Scene_draw_interface;

class VIEWER_EXPORT Viewer : public QGLViewer {

  Q_OBJECT

public:
  Viewer(QWidget * parent, bool antialiasing = false);

  // overload several QGLViewer virtual functions
  void draw();
  void initializeGL();
  void drawWithNames();
  void postSelection(const QPoint&);
  virtual void postDraw();

  void setScene(Scene_draw_interface* scene);
  void setMask(bool b, double ratio=1);
  bool antiAliasing() const { return antialiasing; }

signals:
  void selected(int);
  void pointSelected(QPoint);

public slots:
  void setAntiAliasing(bool b);
  void setTwoSides(bool b);

private:
  void draw_aux(bool with_names);
  void draw_mask();

  Scene_draw_interface* scene;
  bool antialiasing;
  bool twosides;
  bool mask_;
  double ratio_;
}; // end class Viewer

#endif // VIEWER_H
