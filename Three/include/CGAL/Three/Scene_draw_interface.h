#ifndef SCENE_DRAW_INTERFACE_H
#define SCENE_DRAW_INTERFACE_H

class QKeyEvent;
namespace CGAL
{
namespace Three {
  class Viewer_interface;


//! Base class to interact with the scene from the viewer.
class Scene_draw_interface {
public:
  virtual ~Scene_draw_interface(){}
  virtual void initializeGL() = 0;
  virtual void draw() = 0;
  virtual void draw(CGAL::Three::Viewer_interface*) { draw(); };
  virtual void drawWithNames() = 0;
  virtual void drawWithNames(CGAL::Three::Viewer_interface*) { drawWithNames(); }
  virtual bool keyPressEvent(QKeyEvent* e) = 0;
  virtual float get_bbox_length() const = 0;
};
}
}
#endif // SCENE_DRAW_INTERFACE_H;
