#include <CGAL_demo/Scene_item_with_display_list.h>

struct Scene_display_list_item
  : public Scene_item_with_display_list
{
  Scene_display_list_item(GLuint id, GLuint edges_id = 0)
    : Scene_item_with_display_list(id, edges_id) {}

  Scene_item* clone() const { return 0; }
  bool supportsRenderingMode(RenderingMode) const { return false; }
  QString toolTip() const { return tr("Display list item"); }
  void direct_draw() const {}
};
