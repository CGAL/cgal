#ifndef SCENE_TEXTURED_POLYHEDRON_ITEM_H
#define SCENE_TEXTURED_POLYHEDRON_ITEM_H

#include "Scene_textured_polyhedron_item_config.h"
#include "Scene_item_with_display_list.h"
#include "Textured_polyhedron_type_fwd.h"
#include <iostream>
#include "texture.h"

// This class represents a textured polyhedron in the OpenGL scene
class SCENE_TEXTURED_POLYHEDRON_ITEM_EXPORT Scene_textured_polyhedron_item 
  : public Scene_item_with_display_list {
  Q_OBJECT
public:  
  Scene_textured_polyhedron_item();
//   Scene_textured_polyhedron_item(const Scene_textured_polyhedron_item&);
  Scene_textured_polyhedron_item(const Textured_polyhedron& p);
  Scene_textured_polyhedron_item(Textured_polyhedron* const p);
  ~Scene_textured_polyhedron_item();

  Scene_textured_polyhedron_item* clone() const;
  
  // IO
  bool load(std::istream& in);
  bool save(std::ostream& out) const;

  // Function for displaying meta-data of the item
  virtual QString toolTip() const;

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode /* m */) const { return true; }
  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  virtual void direct_draw() const;

  // Get wrapped textured_polyhedron
  Textured_polyhedron*       textured_polyhedron();
  const Textured_polyhedron* textured_polyhedron() const;

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

private:
  Textured_polyhedron* poly;
  Texture texture;

}; // end class Scene_textured_polyhedron_item

#endif // SCENE_TEXTURED_POLYHEDRON_ITEM_H
