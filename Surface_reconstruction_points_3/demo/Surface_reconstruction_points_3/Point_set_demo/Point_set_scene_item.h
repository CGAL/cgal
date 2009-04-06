#ifndef POINT_SET_ITEM_H
#define POINT_SET_ITEM_H

#include "Point_set_scene_item_config.h"
#include "Scene_item_with_display_list.h"
#include "Point_set_demo_types_fwd.h"

#include <iostream>


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

  // OpenGL drawing
  virtual void direct_draw() const;

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
