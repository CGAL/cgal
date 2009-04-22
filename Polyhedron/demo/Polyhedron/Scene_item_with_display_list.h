#ifndef SCENE_ITEM_WITH_DISPLAY_LIST_H
#define SCENE_ITEM_WITH_DISPLAY_LIST_H

#include "Scene_item.h"
#include <CGAL/gl.h>

// This class represents an object in the scene with an OpenGL rendering using display lists
class SCENE_ITEM_EXPORT Scene_item_with_display_list 
  : public Scene_item
{
public:
  Scene_item_with_display_list();
//   Scene_item_with_display_list(const Scene_item_with_display_list&);
  ~Scene_item_with_display_list();

  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  virtual void direct_draw() const = 0;
  // OpenGL drawing using a display list
  virtual void draw() const; 

public slots:
  // Call that once you have finished changing something in the item
  // (either the properties or internal data).
  virtual void changed();

private:
  // display lists
  mutable GLuint display_list;
  mutable bool display_list_built;

}; // end class Scene_item_with_display_list 

#endif // SCENE_ITEM_WITH_DISPLAY_LIST_H
