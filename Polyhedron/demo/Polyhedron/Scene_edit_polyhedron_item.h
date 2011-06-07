#ifndef SCENE_EDIT_POLYHEDRON_ITEM_H
#define SCENE_EDIT_POLYHEDRON_ITEM_H

#include "Scene_edit_polyhedron_item_config.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type_fwd.h"
#include <iostream>

#include <vector>

#include <QColor>

class QMenu;
struct Scene_edit_polyhedron_item_priv;

// This class represents a polyhedron in the OpenGL scene
class SCENE_EDIT_POLYHEDRON_ITEM_EXPORT Scene_edit_polyhedron_item 
  : public Scene_item {
  Q_OBJECT
public:  
  Scene_edit_polyhedron_item(Scene_polyhedron_item* poly_item);
  ~Scene_edit_polyhedron_item();

  Scene_edit_polyhedron_item* clone() const;
  
  // // IO
  // bool load(std::istream& in);
  // bool save(std::ostream& out) const;

  // Function for displaying meta-data of the item
  virtual QString toolTip() const;

  // // Function to override the context menu
  // QMenu* contextMenu();
  
  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode) const { return true; }
  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  void draw() const;

  // Get wrapped polyhedron
  Polyhedron*       polyhedron();
  const Polyhedron* polyhedron() const;

  Scene_polyhedron_item* to_polyhedron_item() const;

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

public slots:
  virtual void changed();
  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z);

protected:
  Scene_edit_polyhedron_item_priv* d;

}; // end class Scene_edit_polyhedron_item

#endif // SCENE_EDIT_POLYHEDRON_ITEM_H
