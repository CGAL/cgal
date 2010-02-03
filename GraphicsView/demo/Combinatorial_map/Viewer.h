#ifndef VIEWER_H
#define VIEWER_H

#include "typedefs.h"
#include <QGLViewer/qglviewer.h>


class Viewer : public QGLViewer {
  Q_OBJECT

  CGAL::Timer timer;
  Scene* scene;
  bool wireframe;
  bool flatShading;
  bool edges;
  bool vertices;

  typedef Map::Dart_handle Dart_handle;

public:
  Viewer(QWidget* parent)
    : QGLViewer(parent), wireframe(false), flatShading(true), edges(true), vertices(true)
  {}

  void setScene(Scene* scene_)
  {
    scene = scene_;
  }

  void clear();

public:
  void draw();
  void Viewer::drawFace(Dart_handle ADart, int AMark);
  void Viewer::drawEdges(Dart_handle ADart);

  virtual void init();
  // void  gl_draw_surface();

  virtual QString helpString() const
  {
    QString text("<h2>M a p   V i e w e r</h2>");
    text += "Use the mouse to move the camera around the object. ";
    text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
    text += "Left and middle buttons pressed together rotate around the camera view direction axis<br><br>";
    text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
    text += "Simply press the function key again to restore it. Several keyFrames define a ";
    text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
    text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
    text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
    text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
    text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
    text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
    text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
    text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
    text += "Press <b>Escape</b> to exit the viewer.";
    return text;
  }
 public slots :

 void sceneChanged();
   
};

#endif
