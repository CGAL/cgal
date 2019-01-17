#ifndef VIEWER_H_
#define VIEWER_H_
#include <QMap>
#include <CGAL/Qt/qglviewer.h>

// forward declarations
class QWidget;
class Scene;
class Viewer : public CGAL::QGLViewer {

  Q_OBJECT

 public:
  explicit Viewer(QWidget *parent);
  Viewer(const Viewer &) = delete;
  Viewer &operator = (const Viewer &) = delete;
  Viewer(const Viewer &&) = delete;
  Viewer &operator = (const Viewer &&) = delete;

  // overload several QGLViewer virtual functions
  void draw();
  void initializeGL();
  void setScene(Scene *pScene);

 protected:
  /*virtual void mousePressEvent(QMouseEvent *e);
  virtual void mouseReleaseEvent(QMouseEvent *e);*/

 private:
  Scene *m_pScene;
  bool m_custom_mouse;
}; // end class Viewer

#endif // VIEWER_H_
