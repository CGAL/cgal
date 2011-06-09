#ifndef SCENE_EDIT_POLYHEDRON_ITEM_H
#define SCENE_EDIT_POLYHEDRON_ITEM_H

#include "Scene_edit_polyhedron_item_config.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include <iostream>

#include <vector>

#include <QColor>
#include <QList>

class QMenu;
struct Scene_edit_polyhedron_item_priv;

// This class represents a polyhedron in the OpenGL scene
class SCENE_EDIT_POLYHEDRON_ITEM_EXPORT Scene_edit_polyhedron_item 
  : public Scene_item {
  Q_OBJECT
public:  
  /// Create an Scene_edit_polyhedron_item from a Scene_polyhedron_item.
  /// The ownership of the polyhedron is moved to the new edit_polyhedron
  /// item.
  Scene_edit_polyhedron_item(Scene_polyhedron_item* poly_item);
  ~Scene_edit_polyhedron_item();

  /// Returns 0, so that one cannot clone an "edit polyhedron" item.
  Scene_edit_polyhedron_item* clone() const;
  
  // // IO
  // bool load(std::istream& in);
  // bool save(std::ostream& out) const;

  // Function for displaying meta-data of the item
  QString toolTip() const;

  // // Function to override the context menu
  // QMenu* contextMenu();
  
  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode) const { return true; }
  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  void draw() const;
 
  bool manipulatable() const { return true; }
  qglviewer::ManipulatedFrame* manipulatedFrame();

  // Get wrapped polyhedron
  Polyhedron*       polyhedron();
  const Polyhedron* polyhedron() const;

  // Functions related to the edition
  Kernel::Point_3 original_position() const;
  Kernel::Point_3 current_position() const;
  Polyhedron::Vertex_handle selected_vertex() const;
  QList<Polyhedron::Vertex_handle> selected_vertices() const;

  /// Returns a Scene_polyhedron_item from the edit polyhedron item, and
  /// transfer the ownership of the polyhedron to it.
  /// The item 'this' must be destroy just after a call to this function.
  Scene_polyhedron_item* to_polyhedron_item() const;

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

public slots:
  void changed();
  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z);
  void setZoneSize(int i);
  void vertex_has_been_selected(void* vertex_handle);

signals:
  void begin_edit();
  void modified();

protected:
  Scene_edit_polyhedron_item_priv* d;

}; // end class Scene_edit_polyhedron_item

#endif // SCENE_EDIT_POLYHEDRON_ITEM_H
