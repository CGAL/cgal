#ifndef _SURFACE_H
#define _SURFACE_H

#include <QObject>
#include <QString>

class Surface : public QObject
{
  Q_OBJECT
public slots:
  virtual void open(const QString& filename) = 0;
  virtual void close() = 0;
  virtual void draw() = 0;
  virtual void get_bbox(float&, float&, float&,
			float&, float&, float&) = 0;
};

#endif // _SURFACE_H
