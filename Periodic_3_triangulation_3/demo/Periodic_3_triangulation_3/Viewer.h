#ifndef VIEWER_H
#define VIEWER_H
#include <QMap>
#include <QGLViewer/qglviewer.h>

class Viewer : public QGLViewer{
    Q_OBJECT

  public:
      Viewer(QWidget* parent);
      virtual ~Viewer();
};
#endif //VIEWER_H
