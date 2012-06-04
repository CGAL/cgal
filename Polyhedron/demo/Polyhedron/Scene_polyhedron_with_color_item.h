#ifndef SCENE_POLYHEDRON_WITH_COLOR_ITEM_H
#define SCENE_POLYHEDRON_WITH_COLOR_ITEM_H

#include "Scene_polyhedron_with_color_item_config.h"
#include "Scene_item_with_display_list.h"
#include "Polyhedron_type_fwd.h"
#include <iostream>

#include <set>
#include <vector>

#include <QColor>

#include <CGAL/Surface_mesh_segmentation.h>
class QMenu;

// This class represents a polyhedron in the OpenGL scene
class SCENE_POLYHEDRON_WITH_COLOR_ITEM_EXPORT Scene_polyhedron_with_color_item 
  : public Scene_item_with_display_list {
  Q_OBJECT
public:  
  Scene_polyhedron_with_color_item();
//   Scene_polyhedron_item(const Scene_polyhedron_item&);
  Scene_polyhedron_with_color_item(const Polyhedron& p);
  Scene_polyhedron_with_color_item(Polyhedron* const p);
  ~Scene_polyhedron_with_color_item();

  Scene_polyhedron_with_color_item* clone() const;
  
  // IO
  bool load(std::istream& in);
  bool save(std::ostream& out) const;

  // Function for displaying meta-data of the item
  virtual QString toolTip() const;

  // Function to override the context menu
  QMenu* contextMenu();
  
  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode) const { return true; }
  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  virtual void direct_draw() const;
  virtual void direct_draw_edges() const;
  virtual void draw() const;
  // Get wrapped polyhedron
  Polyhedron*       polyhedron();
  const Polyhedron* polyhedron() const;

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

  void set_color_vector(std::vector<QColor>& colors);
  
public slots:
  virtual void changed();
  void show_only_feature_edges(bool);
  void enable_facets_picking(bool);
  void set_erase_next_picked_facet(bool);

  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z);

signals:
  void selected_vertex(void*);
  void selected_facet(void*);

private:
  // Initialization
  void init();
  
private:
  Polyhedron* poly;
public:
  CGAL::Surface_mesh_segmentation<Polyhedron>* segmentation;
private:
  typedef Scene_item_with_display_list Base;
  typedef std::vector<QColor> Color_vector;
    
  Color_vector colors_;

  bool show_only_feature_edges_m;
  bool facet_picking_m;
  bool erase_next_picked_facet_m;
}; // end class Scene_polyhedron_item

#endif // SCENE_POLYHEDRON_WITH_COLOR_ITEM_H
