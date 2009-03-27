#ifndef SCENE_POLYHEDRON_ITEM_H
#define SCENE_POLYHEDRON_ITEM_H

#include "Scene_polyhedron_item_config.h"
#include "Scene_item_with_display_list.h"
#include "Polyhedron_type_fwd.h"
#include <iostream>

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
  bool load(std::istream& in);
  bool save(std::ostream& out) const;

  QString toolTip() const;

  void direct_draw() const;
  Polyhedron* polyhedron();

  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

private:
  Polyhedron* poly;
}; // end class Scene_polyhedron_item

#endif // SCENE_POLYHEDRON_ITEM_H
