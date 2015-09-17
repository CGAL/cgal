#ifndef VIEWER_INTERFACE_H
#define VIEWER_INTERFACE_H

#include <QGLViewer/qglviewer.h>
#include <QWidget>
#include <QPoint>
#include <QOpenGLFunctions_2_1>
#include <CGAL/Qt/CreateOpenGLContext.h>
// forward declarations
class QWidget;
class Scene_draw_interface;
class QMouseEvent;
class QKeyEvent;

#include "../Viewer_config.h" // for VIEWER_EXPORT

class VIEWER_EXPORT Viewer_interface : public QGLViewer, public QOpenGLFunctions_2_1 {

  Q_OBJECT

public:
  Viewer_interface(QWidget* parent) : QGLViewer(CGAL::Qt::createOpenGLContext(), parent) {}
  virtual ~Viewer_interface() {}

  virtual void setScene(Scene_draw_interface* scene) = 0;
  virtual bool antiAliasing() const = 0;

  // Those two functions are defined in Viewer.cpp
  static bool readFrame(QString, qglviewer::Frame&);
  static QString dumpFrame(const qglviewer::Frame&);

  virtual bool inFastDrawing() const = 0;
  typedef void (APIENTRYP PFNGLDRAWARRAYSINSTANCEDARBPROC) (GLenum mode, GLint first, GLsizei count, GLsizei primcount);
  typedef void (APIENTRYP PFNGLVERTEXATTRIBDIVISORARBPROC) (GLuint index, GLuint divisor);
  typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTURE2DEXTPROC) (GLuint target, GLuint attachment, GLuint textarget, GLuint texture, GLint level);
  PFNGLDRAWARRAYSINSTANCEDARBPROC glDrawArraysInstanced;
  PFNGLVERTEXATTRIBDIVISORARBPROC glVertexAttribDivisor;
  PFNGLFRAMEBUFFERTEXTURE2DEXTPROC glFramebufferTexture2D;
  bool extension_is_found;
  GLfloat pickMatrix_[16];

Q_SIGNALS:
  void selected(int);
  void requestContextMenu(QPoint global_pos);
  void selectedPoint(double, double, double);
  void selectionRay(double, double, double, double, double, double);

public Q_SLOTS:
  virtual void setAntiAliasing(bool b) = 0;
  virtual void setTwoSides(bool b) = 0;

  virtual void turnCameraBy180Degres() = 0;

  virtual QString dumpCameraCoordinates() = 0;
  virtual bool moveCameraToCoordinates(QString, 
                                       float animation_duration = 0.5f) = 0;

}; // end class Viewer_interface

#endif // VIEWER_INTERFACE_H
