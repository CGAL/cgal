#ifndef _SURFACE_H
#define _SURFACE_H

#include <QObject>
#include <QString>
#include <QGLWidget>

class Viewer;

class Surface : public QObject
{
  Q_OBJECT
protected:
  Surface(QObject* parent) : m_inverse_normals(false) 
  {
    QGLWidget* widget = qFindChild<QGLWidget*>(parent, "viewer");
    if(widget)
      connect(this, SIGNAL(changed()), widget, SLOT(updateGL()));
  }
public slots:
  virtual void open(const QString& filename) = 0;
  virtual void close() = 0;
  virtual void draw() = 0;
  virtual void get_bbox(float&, float&, float&,
			float&, float&, float&) = 0;
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
};

#endif // _SURFACE_H
