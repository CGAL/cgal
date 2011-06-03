#ifndef VIEWER_H
#define VIEWER_H

#include "Viewer_config.h"

#include <QGLViewer/qglviewer.h>
#include <QPoint>

// forward declarations
class QWidget;
class Scene_draw_interface;
class QMouseEvent;
class QKeyEvent;

class Viewer_impl;

class VIEWER_EXPORT Viewer : public QGLViewer {

  Q_OBJECT

public:
  Viewer(QWidget * parent, bool antialiasing = false);
  ~Viewer();

  // overload several QGLViewer virtual functions
  void draw();
  void initializeGL();
  void drawWithNames();
  void postSelection(const QPoint&);

  void setScene(Scene_draw_interface* scene);
  bool antiAliasing() const;

  static bool readFrame(QString, qglviewer::Frame&);
  static QString dumpFrame(const qglviewer::Frame&);

signals:
  void selected(int);
  void requestContextMenu(QPoint global_pos);
  void selectedPoint(double, double, double);
  void selectionRay(double, double, double, double, double, double);

public slots:
  void setAntiAliasing(bool b);
  void setTwoSides(bool b);

  void turnCameraBy180Degres();

  QString dumpCameraCoordinates();
  bool moveCameraToCoordinates(QString, 
                               float animation_duration = 0.5f);

protected:
  void mousePressEvent(QMouseEvent*);
  void keyPressEvent(QKeyEvent*);

protected:
  Viewer_impl* d;
}; // end class Viewer

#endif // VIEWER_H
