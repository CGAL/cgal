#ifndef _VIEWER_H
#define _VIEWER_H

#include <QGLViewer/qglviewer.h>

class Viewer : public QGLViewer
{
public:
   Viewer(QWidget* parent) : QGLViewer(parent) {};
protected :
  virtual void init();
  virtual QString helpString() const;
};

#endif // _VIEWER_H
