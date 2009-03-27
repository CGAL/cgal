#include "Scene_item_with_display_list.h"
#include <iostream>

Scene_item_with_display_list::Scene_item_with_display_list()
  : display_list(0),
    display_list_built(false)
{}

// Scene_item_with_display_list::
// Scene_item_with_display_list(const Scene_item_with_display_list& item)
//   : Scene_item(item),
//     display_list(0),
//     display_list_built(false)
// {}

Scene_item_with_display_list::~Scene_item_with_display_list()
{
  if(display_list_built && display_list != 0) {
    ::glDeleteLists(display_list,1);
  }
}

void Scene_item_with_display_list::draw() const
{
  if(!display_list_built)
  {
    if(display_list == 0) {
      display_list = ::glGenLists(1);
      if(display_list == 0)
      {
        std::cerr << "Unable to create display list" << std::endl;
        return;
      }
    }
    // draw the mesh in a display list
    ::glNewList(display_list,GL_COMPILE_AND_EXECUTE);
    direct_draw();
    ::glEndList();
    display_list_built = true;
  }
  else {
    ::glCallList(display_list);
  }
}

void Scene_item_with_display_list::changed()
{
  display_list_built = false;
}
