#ifndef SCENE_DRAW_INTERFACE_H
#define SCENE_DRAW_INTERFACE_H

class QKeyEvent;
class Viewer_interface;

class Scene_draw_interface {
public:
  virtual ~Scene_draw_interface(){}
  virtual void initializeGL() = 0;
  virtual void draw() = 0;
  virtual void draw(Viewer_interface*) { draw(); };
  virtual void drawWithNames() = 0;
  virtual void drawWithNames(Viewer_interface*) { drawWithNames(); }
  virtual bool keyPressEvent(QKeyEvent* e) = 0;
};

#endif // SCENE_DRAW_INTERFACE_H;
