#include <CGAL_demo/Scene_item_with_display_list.h>
#include <iostream>

Scene_item_with_display_list::Scene_item_with_display_list()
{
  for(int i = 0; i < NB_OF_DISPLAY_LISTS; ++i) {
    display_list[i] = 0;
    display_list_built[i] = false;
  }
}

// Scene_item_with_display_list::
// Scene_item_with_display_list(const Scene_item_with_display_list& item)
//   : Scene_item(item),
//     display_list(0),
//     display_list_built(false)
// {}

Scene_item_with_display_list::~Scene_item_with_display_list()
{
  for(int i = 0; i < NB_OF_DISPLAY_LISTS; ++i) {
    if(display_list_built[i] && display_list[i] != 0) {
      ::glDeleteLists(display_list[i],1);
    }
  }
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list

void Scene_item_with_display_list::draw(QGLViewer *viewer)  {
  draw(viewer, DRAW);
}

void Scene_item_with_display_list::draw_edges(QGLViewer *viewer)  {
  draw(viewer, DRAW_EDGES);
}

void Scene_item_with_display_list::draw(QGLViewer *viewer, int i)
{
  if(!display_list_built[i])
  {
    if(display_list[i] == 0) {
      display_list[i] = ::glGenLists(1);
      if(display_list[i] == 0)
      {
        std::cerr << "Unable to create display list" << std::endl;
        return;
      }
    }
    // draw the item in a display list
    ::glNewList(display_list[i],GL_COMPILE_AND_EXECUTE);
    if(i == 0) {
      direct_draw(viewer);
    }
    else {
      direct_draw_edges(viewer);
    }
    ::glEndList();
    display_list_built[i] = true;
  }
  else {
    // draw using the display list
    ::glCallList(display_list[i]);
  }
}

void Scene_item_with_display_list::changed()
{
  for(int i = 0; i < NB_OF_DISPLAY_LISTS; ++i) {
    display_list_built[i] = false;
  }


}
