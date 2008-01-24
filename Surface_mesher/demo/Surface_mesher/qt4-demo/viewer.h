#ifndef _VIEWER_H
#define _VIEWER_H

#include <QGLViewer/qglviewer.h>

class Surface;

class Viewer : public QGLViewer
{
  Q_OBJECT
public:
  Viewer(QWidget* parent) : QGLViewer(parent), parent(parent), surface(0) {};

  void set_surface(Surface*);

public slots:
  void interpolateToFitBoundingBox(double, double, double, double, double, double);

protected :
  virtual void init();
  virtual void draw();
  virtual void drawWithNames();
  virtual void postSelection(const QPoint&);
  virtual QString helpString() const;

  QWidget* parent;
  Surface* surface;
};

#endif // _VIEWER_H
