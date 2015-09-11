#ifndef SCENE_POLYHEDRON_ITEM_DECORATOR_H
#define SCENE_POLYHEDRON_ITEM_DECORATOR_H
#include "Scene_polyhedron_item_decorator_config.h"
#include "Scene_polyhedron_item.h"

// This class is a decorator for Scene_polyhedron_item yet it does not inherit it but Scene_item
class SCENE_POLYHEDRON_ITEM_DECORATOR_EXPORT Scene_polyhedron_item_decorator 
  : public Scene_item {
  Q_OBJECT
public:  
  /// Create an Scene_polyhedron_item_decorator from a Scene_polyhedron_item.

  Scene_polyhedron_item_decorator(Scene_polyhedron_item* poly_item, bool delete_item = true);
  ~Scene_polyhedron_item_decorator();

  /// Returns 0, so that one cannot clone decorator
  Scene_polyhedron_item_decorator* clone() const;
  
  // // IO
  // bool load(std::istream& in);
  // bool save(std::ostream& out) const;

  // Function for displaying meta-data of the item
  QString toolTip() const;

  // // Function to override the context menu
  // QMenu* contextMenu();
  
  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const { return (m!=PointsPlusNormals && m!=Splatting); }
  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  // dispatch to poly_item direct_draw and direct_draw_edges
  void draw() const;
  void draw_edges() const;

  // Get wrapped polyhedron
  Polyhedron*       polyhedron();
  const Polyhedron* polyhedron() const;

  Scene_polyhedron_item* polyhedron_item() const;
  void                   set_polyhedron_item(Scene_polyhedron_item* poly_item);

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

  bool delete_item() { return delete_poly_item; }
  void set_delete_item(bool delete_item) { delete_poly_item = delete_item; }

public Q_SLOTS:
  void invalidate_buffers();
  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z);

protected:
  Scene_polyhedron_item* poly_item;
  bool delete_poly_item;
}; // end class Scene_polyhedron_item_decorator

#endif // SCENE_POLYHEDRON_ITEM_DECORATOR_H
