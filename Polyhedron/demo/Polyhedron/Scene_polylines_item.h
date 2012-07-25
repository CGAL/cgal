#ifndef SCENE_POLYLINES_ITEM_H
#define SCENE_POLYLINES_ITEM_H

#include "Scene_polylines_item_config.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "Scene_item.h"

#include <QString>

#include <list>
#include <vector>

class Scene_polylines_item_private;

class SCENE_POLYLINES_ITEM_EXPORT Scene_polylines_item : public Scene_item
{
  Q_OBJECT
public:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_3 Point_3;
  typedef std::vector<Point_3> Polyline;
  typedef std::list<Polyline> Polylines_container;

  typedef K::Iso_cuboid_3 Iso_cuboid_3;

  Scene_polylines_item();
  virtual ~Scene_polylines_item();

  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

  Scene_polylines_item* clone() const;

  QString toolTip() const;

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const;

  QMenu* contextMenu();
  
  // Flat/Gouraud OpenGL drawing
  void draw() const;

  // Wireframe OpenGL drawing
  void draw_edges() const;

  void draw_points() const;

public slots:
  void change_corner_radii(double);
  void change_corner_radii();
  void split_at_sharp_angles();

  void merge(Scene_polylines_item*);

public:
  Polylines_container polylines;

  // http://en.wikipedia.org/wiki/D-pointer
  Scene_polylines_item_private* d;

}; // end class Scene_polylines_item

#endif
