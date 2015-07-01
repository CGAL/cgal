#ifndef _VIEWER_H
#define _VIEWER_H

#include <QGLViewer/qglviewer.h>

class Surface;

class Viewer : public QGLViewer
{
  Q_OBJECT
public:
  Viewer(QWidget* parent);

  void set_surface(Surface*);

public Q_SLOTS:
  void interpolateToFitBoundingBox(double, double, double, double, double, double);

protected :
  virtual void init();
  virtual void draw();
  virtual void drawWithNames();
  virtual void postSelection(const QPoint&);
  virtual QString helpString() const;

  Surface* surface;
};

#endif // _VIEWER_H
