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
  Surface(QObject* parent) : m_inverse_normals(false) 
  {
    viewer = qFindChild<Viewer*>(parent, "viewer");
    if(viewer)
      connect(this, SIGNAL(changed()), viewer, SLOT(updateGL()));
    viewer->set_surface(this);
  }
public slots:
  virtual void open(const QString& filename) = 0;
  virtual void close() = 0;
  virtual void draw() = 0;
  virtual void get_bbox(float&, float&, float&,
			float&, float&, float&) = 0;
  virtual void drawWithNames() {};
  virtual void postSelection(const QPoint&) {};
signals:
  void changed();

public slots:
  void set_inverse_normals(const bool b) {
    m_inverse_normals = b;
    emit(changed());
  }

public:
  bool inverse_normals() const {
    return m_inverse_normals;
  }

protected:
  bool m_inverse_normals;
  Viewer* viewer;
};

#endif // _SURFACE_H
