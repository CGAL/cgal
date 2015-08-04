#ifndef VIEWER_H
#define VIEWER_H

#include "Viewer_config.h"
#include <CGAL_demo/Viewer_interface.h>

#include <QGLViewer/qglviewer.h>
#include <QPoint>

// forward declarations
class QWidget;
class Scene_draw_interface;
class QMouseEvent;
class QKeyEvent;

class Viewer_impl;

class VIEWER_EXPORT Viewer : public Viewer_interface {

  Q_OBJECT

public:
  Viewer(QWidget * parent, bool antialiasing = false);
  ~Viewer();

  // overload several QGLViewer virtual functions
  void draw();
  void fastDraw();
  void initializeGL();
  void drawWithNames();
  void postSelection(const QPoint&);
  void beginSelection(const QPoint &point);
  void endSelection(const QPoint &point);
  void setScene(Scene_draw_interface* scene);
  bool antiAliasing() const;

  bool inFastDrawing() const;

public Q_SLOTS:
  void setAntiAliasing(bool b);
  void setTwoSides(bool b);

  void turnCameraBy180Degres();

  QString dumpCameraCoordinates();
  bool moveCameraToCoordinates(QString, 
                               float animation_duration = 0.5f);

protected:
  void mousePressEvent(QMouseEvent*);
  void keyPressEvent(QKeyEvent*);
  void pickMatrix(GLdouble x, GLdouble y, GLdouble width, GLdouble height,
                  GLint viewport[4]);

protected:
  Viewer_impl* d;
}; // end class Viewer

#endif // VIEWER_H
