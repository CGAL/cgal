#ifndef SCENE_COMBINATORIAL_MAP_ITEM_H
#define SCENE_COMBINATORIAL_MAP_ITEM_H


//=========
#include <CGAL/internal/corefinement/Combinatorial_map_for_corefinement.h>
#include "Scene_combinatorial_map_item_config.h"
#include "Scene_item_with_display_list.h"
#include <iostream>

#include "Polyhedron_type.h"

typedef CGAL::internal_IOP::Item_with_points_and_volume_info<Kernel,Polyhedron> Items;
typedef CGAL::Combinatorial_map<3,Items> Combinatorial_map_3;
//=========

class QMenu;
class QAction;
class Scene_interface;
class Scene_polyhedron_item;

class SCENE_COMBINATORIAL_MAP_ITEM_EXPORT Scene_combinatorial_map_item
  : public Scene_item_with_display_list 
{
  Q_OBJECT
public:  
  Scene_combinatorial_map_item(Scene_interface*,void* ad_A=NULL);
  ~Scene_combinatorial_map_item();

  Scene_combinatorial_map_item* clone() const;
  // Function to override the context menu
  QMenu* contextMenu();

//  bool load(std::istream& in);
//  void load(Scene_polyhedron_item*);
//  bool save(std::ostream& out) const;

  QString toolTip() const;

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const { return (m != Gouraud && m!=PointsPlusNormals); } // CHECK THIS!
  //Event handling
  virtual bool keyPressEvent(QKeyEvent*);
  // OpenGL drawing in a display list
  void direct_draw() const;
  void direct_draw_edges() const;
  void draw_points() const;

  bool isFinite() const { return true; }
  bool is_from_corefinement() const {return address_of_A!=NULL;}
  bool isEmpty() const;
  Bbox bbox() const;
  
  const Combinatorial_map_3& combinatorial_map() const
  {
    return *m_combinatorial_map;
  }

  Combinatorial_map_3& combinatorial_map()
  {
    return *m_combinatorial_map;
  }
  
  Combinatorial_map_3* m_combinatorial_map;
  
private:
  Kernel::Vector_3 compute_face_normal(Combinatorial_map_3::Dart_const_handle adart) const;
  Scene_interface* last_known_scene;
  std::size_t volume_to_display;
  QAction* exportSelectedVolume;
  void* address_of_A;
  template <class Predicate> void export_as_polyhedron(Predicate,const QString&) const;

public slots:
  void set_next_volume();
  void export_current_volume_as_polyhedron() const;
  void export_union_as_polyhedron() const;
  void export_intersection_as_polyhedron() const;
  void export_A_minus_B_as_polyhedron() const;
  void export_B_minus_A_as_polyhedron() const;
}; // end class Scene_combinatorial_map_item

#endif // SCENE_COMBINATORIAL_MAP_ITEM_H
