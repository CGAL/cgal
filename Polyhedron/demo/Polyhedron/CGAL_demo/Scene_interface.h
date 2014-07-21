#ifndef SCENE_INTERFACE_H
#define SCENE_INTERFACE_H

#include <QString>
#include <QColor>
#include <QList>
#include <algorithm>
#include <cmath>

class Scene_item;

// OpenGL rendering mode
enum RenderingMode { Points = 0,
                     PointsPlusNormals,
                     Splatting,
                     Wireframe, 
                     Flat,
                     FlatPlusEdges,
                     Gouraud,
                     LastRenderingMode = Gouraud,
                     NumberOfRenderingMode = LastRenderingMode+1 };

// Interface of Scene class exported to plugins
class Scene_interface {
public:
  struct Bbox {
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;
    Bbox(const double _xmin,const double _ymin,const double _zmin,
         const double _xmax,const double _ymax,const double _zmax)
      : xmin(_xmin), ymin(_ymin), zmin(_zmin),
        xmax(_xmax), ymax(_ymax), zmax(_zmax)
    {
    }
    Bbox()
      : xmin(0.0), ymin(0.0), zmin(0.0),
        xmax(1.0), ymax(1.0), zmax(1.0)
    {
    }

    Bbox operator+(const Bbox& b) const {
      return Bbox((std::min)(xmin, b.xmin),
                  (std::min)(ymin, b.ymin),
                  (std::min)(zmin, b.zmin),
                  (std::max)(xmax, b.xmax),
                  (std::max)(ymax, b.ymax),
                  (std::max)(zmax, b.zmax));
    }

    bool operator==(const Bbox&b) const{
      return
        xmin==b.xmin && xmax==b.xmax &&
        ymin==b.ymin && ymax==b.ymax &&
        zmin==b.zmin && zmax==b.zmax;
    }

    bool operator!=(const Bbox& b) const{
      return !(*this == b);
    }
    
    double width() const { return xmax-xmin; }
    double height() const { return ymax-ymin; }
    double depth() const { return zmax-zmin; }
    
    double diagonal_length() const
    {
      return std::sqrt(width()*width() + height()*height() + depth()*depth());
    }

  }; // struct BBox (ad hoc class, does not depend on CGAL kernels

  typedef int Item_id;

  virtual ~Scene_interface() {};

  virtual Item_id addItem(Scene_item* item) = 0;
  virtual Scene_item* replaceItem(Item_id, Scene_item*, bool emit_item_about_to_be_destroyed = false) = 0;

  virtual Item_id erase(Item_id) = 0;
  // Returns the index of the item just before the one that is erased,
  // or just after. Returns -1 if the list is empty.

  virtual Item_id duplicate(Item_id) = 0;
  // Returns the index of the new item
  // If no new item has been created (because the item type is note
  // clonable), returns -1.

  // Accessors (getters)
  virtual int numberOfEntries() const = 0;
  virtual Scene_item* item(Item_id) const = 0;
  virtual Item_id item_id(Scene_item*) const = 0;
  virtual Item_id mainSelectionIndex() const = 0;
  virtual QList<Item_id> selectionIndices() const = 0;
  virtual Item_id selectionAindex() const = 0;
  virtual Item_id selectionBindex() const = 0;

  // Get scene bounding box
  virtual Bbox bbox() const = 0;
  virtual double len_diagonal() const = 0;

public:
  // Notify the scene that an item was modified
  virtual void itemChanged(Item_id i) = 0; 
  virtual void itemChanged(Scene_item*) = 0;

  // Select an item
  virtual void setSelectedItem(Item_id) = 0;
  
}; // end interface Scene_interface


#endif // SCENE_INTERFACE_H
