#ifndef SCENE_POLYHEDRON_ITEM_H
#define SCENE_POLYHEDRON_ITEM_H

#include "Scene_polyhedron_item_config.h"
#include "Scene_item_with_display_list.h"
#include "Polyhedron_type_fwd.h"
#include <iostream>

// This class represents a polyhedron in the OpenGL scene
class SCENE_POLYHEDRON_ITEM_EXPORT Scene_polyhedron_item
  : public Scene_item_with_display_list {
  Q_OBJECT
public:
  Scene_polyhedron_item();
//   Scene_polyhedron_item(const Scene_polyhedron_item&);
  Scene_polyhedron_item(const Polyhedron& p);
  Scene_polyhedron_item(Polyhedron* const p);
  ~Scene_polyhedron_item();

  Scene_polyhedron_item* clone() const;

  // IO
  bool load(std::istream& in);
  bool save(std::ostream& out) const;

  // Function for displaying meta-data of the item
  virtual QString toolTip() const;

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const;
  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  virtual void direct_draw() const;

  // Gets wrapped polyhedron
  Polyhedron*       polyhedron();
  const Polyhedron* polyhedron() const;

  // Gets dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

private:
  Polyhedron* poly;

}; // end class Scene_polyhedron_item

#endif // SCENE_POLYHEDRON_ITEM_H
