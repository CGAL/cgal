#ifndef _SURFACE_H
#define _SURFACE_H

#include <QObject>
#include <QString>
#include <QGLWidget>

#include "viewer.h"

class Surface : public QObject
{
  Q_OBJECT
protected:
  Surface(QObject* parent)
  {
    viewer = qFindChild<Viewer*>(parent, "viewer");
    if(viewer)
      connect(this, SIGNAL(changed()), viewer, SLOT(updateGL()));
    viewer->set_surface(this);
  }
public slots:
  virtual bool open(const QString& filename) = 0;
  virtual void close() = 0;
  virtual void draw() = 0;
  virtual void get_bbox(float&, float&, float&,
			float&, float&, float&) = 0;
  virtual void drawWithNames() {};
  virtual void postSelection(const QPoint&) {};
signals:
  void changed();

protected:
  Viewer* viewer;
};

#endif // _SURFACE_H
