#ifndef SCENE_POLYGON_SOUP_ITEM_H
#define SCENE_POLYGON_SOUP_ITEM_H

#include "Scene_polygon_soup_item_config.h"
#include "Scene_item_with_display_list.h"
#include <iostream>

#include "Polyhedron_type_fwd.h"

struct Polygon_soup;
class Scene_polyhedron_item;

class SCENE_POLYGON_SOUP_EXPORT Scene_polygon_soup_item 
  : public Scene_item_with_display_list 
{
  Q_OBJECT
public:  
  Scene_polygon_soup_item();
  ~Scene_polygon_soup_item();

  Scene_polygon_soup_item* clone() const;
  bool load(std::istream& in);
  void load(Scene_polyhedron_item*);
  bool save(std::ostream& out) const;

  QString toolTip() const;

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const { return m != Gouraud; } // CHECK THIS!
  // OpenGL drawing in a display list
  void direct_draw() const;
  void draw_points() const;

  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

  void new_vertex(const double&, const double&, const double&);
  void new_triangle(const std::size_t, const std::size_t, const std::size_t);

public slots:
  void shuffle_orientations();
  bool orient();
  bool exportAsPolyhedron(Polyhedron*);
  void inside_out();

  void setDisplayNonManifoldEdges(const bool);
  bool displayNonManifoldEdges() const;

private:
  Polygon_soup* soup;
  bool oriented;
}; // end class Scene_polygon_soup_item

#endif // SCENE_POLYGON_SOUP_ITEM_H
