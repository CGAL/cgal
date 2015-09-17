#ifndef VIEWER_H
#define VIEWER_H

#include "Viewer_config.h"
#include <CGAL_demo/Viewer_interface.h>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLShaderProgram>


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
  void drawVisualHints();
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
  struct AxisData
  {
      std::vector<float> *vertices;
      std::vector<float> *normals;
      std::vector<float> *colors;
  };
  QOpenGLBuffer buffers[3];
  QOpenGLVertexArrayObject vao[1];
  QOpenGLShaderProgram rendering_program;
  std::vector<float> v_Axis;
  std::vector<float> n_Axis;
  std::vector<float> c_Axis;
  bool axis_are_displayed;
  void mousePressEvent(QMouseEvent*);
  void keyPressEvent(QKeyEvent*);
  void pickMatrix(GLdouble x, GLdouble y, GLdouble width, GLdouble height,
                  GLint viewport[4]);
  void makeArrow(float R, int prec, qglviewer::Vec from, qglviewer::Vec to, qglviewer::Vec color, AxisData &data);
  void resizeGL(int w, int h);


protected:
  Viewer_impl* d;
}; // end class Viewer

#endif // VIEWER_H
