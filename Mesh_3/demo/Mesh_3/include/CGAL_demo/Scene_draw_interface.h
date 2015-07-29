#ifndef SCENE_DRAW_INTERFACE_H
#define SCENE_DRAW_INTERFACE_H

class Viewer;

class Scene_draw_interface {
public:
  virtual ~Scene_draw_interface(){}
  virtual void initializeGL() = 0;
  virtual void draw(Viewer*) = 0;
  virtual void drawWithNames(Viewer*) = 0;
};

#endif // SCENE_DRAW_INTERFACE_H;
