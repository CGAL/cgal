#ifndef SCENE_ITEM_WITH_DISPLAY_LIST_H
#define SCENE_ITEM_WITH_DISPLAY_LIST_H

#include "Scene_item.h"
#include <Scene_interface.h>
#include <CGAL/gl.h>

class Viewer_interface;

// This class represents an object in the scene with an OpenGL rendering using display lists
class SCENE_ITEM_EXPORT Scene_item_with_display_list 
  : public Scene_item
{
public:
  Scene_item_with_display_list();
  Scene_item_with_display_list(GLuint list_id, GLuint edges_list_id = 0);
//   Scene_item_with_display_list(const Scene_item_with_display_list&);
  ~Scene_item_with_display_list();

  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  virtual void direct_draw() const = 0;
  virtual void direct_draw(Viewer_interface*) const { direct_draw(); }
  virtual void direct_draw_edges() const { draw(); };
  virtual void direct_draw_edges(Viewer_interface*) const { direct_draw_edges(); }
  // OpenGL drawing using a display list
  virtual void draw() const;
  virtual void draw(Viewer_interface*) const { draw(); }
  virtual void draw_edges() const;
  virtual void draw_edges(Viewer_interface*) const { draw_edges(); }

public slots:
  // Call that once you have finished changing something in the item
  // (either the properties or internal data).
  virtual void changed();

protected:
  enum { DRAW = 0, DRAW_EDGES = 1, NB_OF_DISPLAY_LISTS = 2};
  void draw(int) const;
  // display lists
  mutable GLuint display_list[NB_OF_DISPLAY_LISTS];
  mutable bool display_list_built[NB_OF_DISPLAY_LISTS];

}; // end class Scene_item_with_display_list 

#endif // SCENE_ITEM_WITH_DISPLAY_LIST_H
