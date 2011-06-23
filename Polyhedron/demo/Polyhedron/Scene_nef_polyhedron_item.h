#ifndef SCENE_NEF_POLYHEDRON_ITEM_H
#define SCENE_NEF_POLYHEDRON_ITEM_H

#include "Scene_nef_polyhedron_item_config.h"
#include "Scene_item_with_display_list.h"
#include "Nef_type_fwd.h"
#include <iostream>

class Scene_polyhedron_item;

class SCENE_NEF_POLYHEDRON_ITEM_EXPORT Scene_nef_polyhedron_item
 : public Scene_item_with_display_list 
{
  Q_OBJECT
public:  
  Scene_nef_polyhedron_item();
//   Scene_nef_polyhedron_item(const Scene_nef_polyhedron_item&);
  Scene_nef_polyhedron_item(const Nef_polyhedron& p);
  Scene_nef_polyhedron_item(Nef_polyhedron* const p);
  ~Scene_nef_polyhedron_item();

  Scene_nef_polyhedron_item* clone() const;
  bool load_from_off(std::istream& in);
  bool load(std::istream& in);
  bool save(std::ostream& in) const;

  QFont font() const;
  QString toolTip() const;

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const { return m != Gouraud; } // CHECK THIS!
  // OpenGL drawing in a display list
  void direct_draw() const;
  // Wireframe OpenGL drawing
  void draw_edges() const;

  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

  Nef_polyhedron* nef_polyhedron();
  const Nef_polyhedron* nef_polyhedron() const;

  bool is_simple() const;

  // conversion operations
  static Scene_nef_polyhedron_item* from_polyhedron(Scene_polyhedron_item*);
  Scene_polyhedron_item* convert_to_polyhedron() const;

  // Nef boolean operations
  Scene_nef_polyhedron_item&
  operator+=(const Scene_nef_polyhedron_item&); // union

  Scene_nef_polyhedron_item&
  operator*=(const Scene_nef_polyhedron_item&); // intersection

  Scene_nef_polyhedron_item&
  operator-=(const Scene_nef_polyhedron_item&); // difference

  static Scene_nef_polyhedron_item*
  sum(const Scene_nef_polyhedron_item&, 
      const Scene_nef_polyhedron_item&);

  void convex_decomposition(std::list< Scene_polyhedron_item*>&);
  
private:
  Nef_polyhedron* nef_poly;
}; // end class Scene_nef_polyhedron_item

#endif // SCENE_NEF_POLYHEDRON_ITEM_H
