#ifndef POINT_SET_ITEM_H
#define POINT_SET_ITEM_H

#include "Point_set_scene_item_config.h"
#include "Polyhedron_type_fwd.h"
#include "Kernel_type.h"
#include "Point_set_3.h"
#include "Scene_item_with_display_list.h"

#include <iostream>


// point set
typedef Point_set_3<Kernel> Point_set;
typedef Point_set::UI_point UI_point; // type of points in Point_set_3


// This class represents a point set in the OpenGL scene
class POINT_SET_ITEM_EXPORT Point_set_scene_item
  : public Scene_item_with_display_list
{
  Q_OBJECT

public:
  Point_set_scene_item();
  Point_set_scene_item(const Point_set_scene_item& toCopy);
  Point_set_scene_item(const Polyhedron& p);
  ~Point_set_scene_item();
  Point_set_scene_item* clone() const;

  // IO
  bool read_off_point_set(std::istream& in);
  bool write_off_point_set(std::ostream& out) const;
  bool read_xyz_point_set(std::istream& in);
  bool write_xyz_point_set(std::ostream& out) const;

  // Function for displaying meta-data of the item
  virtual QString toolTip() const;

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const;
  // Points OpenGL drawing in a display list
  virtual void direct_draw() const;
  // Normals OpenGL drawing
  virtual void draw_normals() const;
  // Draws oriented points with radius
  virtual void draw_splats() const;

  // Get wrapped point set
  Point_set*       point_set();
  const Point_set* point_set() const;

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

// Data
private:
  Point_set* m_points;

}; // end class Point_set_scene_item


#endif // POINT_SET_ITEM_H
