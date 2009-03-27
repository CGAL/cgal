#ifndef SCENE_ITEM_WITH_DISPLAY_LIST_H
#define SCENE_ITEM_WITH_DISPLAY_LIST_H

#include "Scene_item.h"
#include <CGAL/gl.h>

class SCENE_ITEM_EXPORT Scene_item_with_display_list 
  : public Scene_item
{
public:
  Scene_item_with_display_list();
//   Scene_item_with_display_list(const Scene_item_with_display_list&);

  ~Scene_item_with_display_list();

  virtual void direct_draw() const = 0;
  void draw() const; 
  void changed();

private:
  // display lists
  mutable GLuint display_list;
  mutable bool display_list_built;
}; // end class Scene_item

#endif // SCENE_ITEM_WITH_DISPLAY_LIST_H
