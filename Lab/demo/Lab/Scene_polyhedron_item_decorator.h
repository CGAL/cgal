#ifndef SCENE_POLYHEDRON_ITEM_DECORATOR_H
#define SCENE_POLYHEDRON_ITEM_DECORATOR_H
#include "Scene_polyhedron_item_decorator_config.h"

#include "Scene_surface_mesh_item.h"
typedef Scene_surface_mesh_item Scene_face_graph_item;

typedef Scene_face_graph_item::Face_graph Face_graph;

// This class is a decorator for Scene_surface_mesh_item yet it does not inherit it but
// Scene_item_rendering_helper
class SCENE_POLYHEDRON_ITEM_DECORATOR_EXPORT Scene_polyhedron_item_decorator
  : public CGAL::Three::Scene_item_rendering_helper {
  Q_OBJECT
public:
  /// Create an Scene_polyhedron_item_decorator from a Scene_polyhedron_item.

  Scene_polyhedron_item_decorator(Scene_face_graph_item* poly_item, bool delete_item = true);
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
  bool supportsRenderingMode(RenderingMode m) const { return (m!=PointsPlusNormals ); }

  // Get wrapped polyhedron
  Face_graph*       polyhedron();
  const Face_graph* polyhedron() const;

  Scene_face_graph_item* polyhedron_item() const;
  void                   set_polyhedron_item(Scene_face_graph_item* poly_item);

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  void compute_bbox() const;

  bool delete_item() { return delete_poly_item; }
  void set_delete_item(bool delete_item) { delete_poly_item = delete_item; }

public Q_SLOTS:
  void invalidateOpenGLBuffers();
  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z);

protected:
  Scene_face_graph_item* poly_item;
  bool delete_poly_item;
}; // end class Scene_polyhedron_item_decorator

#endif // SCENE_POLYHEDRON_ITEM_DECORATOR_H
