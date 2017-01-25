// Copyright (c) 2012-2015  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent RINEAU

//! \file Scene_interface.h
#ifndef SCENE_INTERFACE_H
#define SCENE_INTERFACE_H

#include <CGAL/license/Three.h>


#include <QString>
#include <QColor>
#include <QList>
#include <algorithm>
#include <cmath>
#include <CGAL/Bbox_3.h>
namespace CGAL{namespace Three{
class Scene_item;
class Scene_group_item;
    }}
// OpenGL rendering mode
/*!
 * Contains the class Scene_interface and the enum RenderingMode.
 */
/*!
 * \brief The RenderingMode enum
 * Can be Points, PointsPlusNormals, Wireframe, Flat, FlatPlusEdges, or Gouraud.
 * - Points, PointsPlusNormals, and Wireframe have no light model.
 * - Flat and FlatPlusEdges use a basic light model with one normal per facet.
 * - Gouraud use the same light model but with one normal per vertex.
 */
enum RenderingMode { Points = 0,
                     PointsPlusNormals,
                     Splatting,
                     Wireframe, 
                     Flat,
                     FlatPlusEdges,
                     Gouraud,
                     ShadedPoints,
                     NumberOfRenderingMode};


namespace CGAL {
namespace Three{
/*!
 * This is the class given to the plugins to interact with the scene.
 * */
class Scene_interface {
public:
  typedef CGAL::Bbox_3 Bbox;
    //!Contains the coordinates for the two extremities of the diagonal of a bounding box.
    //! Those points are respectively the left bottom front vertex and the right upper back vertex.
 /* struct Bbox {
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;
    Bbox(const double _xmin,const double _ymin,const double _zmin,
         const double _xmax,const double _ymax,const double _zmax)
      : xmin(_xmin), ymin(_ymin), zmin(_zmin),
        xmax(_xmax), ymax(_ymax), zmax(_zmax)
    {
    }

    //!The default constructor of a bounding box.
    //! The two vertices are (0.0,0.0,0.0) and (1.0,1.0,1.0).
    Bbox()
      : xmin(0.0), ymin(0.0), zmin(0.0),
        xmax(1.0), ymax(1.0), zmax(1.0)
    {
    }
    //!Sums the coordinates of each vertex one by one.
    Bbox operator+(const Bbox& b) const {
      return Bbox((std::min)(xmin, b.xmin),
                  (std::min)(ymin, b.ymin),
                  (std::min)(zmin, b.zmin),
                  (std::max)(xmax, b.xmax),
                  (std::max)(ymax, b.ymax),
                  (std::max)(zmax, b.zmax));
    }

    //!@returns true if all the coordinates of the two bounding boxes are the same.
    bool operator==(const Bbox&b) const{
      return
        xmin==b.xmin && xmax==b.xmax &&
        ymin==b.ymin && ymax==b.ymax &&
        zmin==b.zmin && zmax==b.zmax;
    }
    //!@returns the opposite of ==.
    //!@see #operator==
    bool operator!=(const Bbox& b) const{
      return !(*this == b);
    }
    //!@returns the difference between xmax and xmin.
    double width() const { return xmax-xmin; }
    //!@returns the difference between ymax and ymin.
    double height() const { return ymax-ymin; }
    //!@returns the difference between zmax and zmin.
    double depth() const { return zmax-zmin; }
    //!@returns the length of the diagonal of the bounding box.
    double diagonal_length() const
    {
      return CGAL::sqrt((xmax()-xmin())*(xmax()-xmin()) + (ymax()-ymin())*(ymax()-ymin()) + (zmax()-zmin())*(zmax()-zmin()));
    }

  }; // struct BBox (ad hoc class, does not depend on CGAL kernels
*/
  //!Integer used as the index of a Scene_item.
  typedef int Item_id;
  //!Virtual destructor
  virtual ~Scene_interface() {};
  //!Adds an item to the list of items.
  //!@returns the index of the new item.
  virtual Item_id addItem(CGAL::Three::Scene_item* item) = 0;
  //!Replace an item by a new one in the scene. The first is deleted and gives its index to the latter.
  //!If emit_item_about_to_be_destroyed is true, emits
  //!an itemAboutToBeDestroyed signal.
  //!@returns a pointer to the old item.
  virtual Scene_item* replaceItem(Item_id, CGAL::Three::Scene_item*, bool emit_item_about_to_be_destroyed = false) = 0;
  //!Moves item to the targeted group.
  virtual void changeGroup(CGAL::Three::Scene_item* item, CGAL::Three::Scene_group_item* target_group) = 0;

  /*! Erases an item in the list.
   * @returns the index of the item just before the one that is erased,
   * or just after.
   * @returns -1 if the list is empty.*/
  virtual Item_id erase(Item_id) = 0;

  /*! Duplicates a scene item.
   * @returns the index of the new item (-1 on error).
   */
    virtual Item_id duplicate(Item_id) = 0;

  // Accessors (getters)
  //!
  //!@returns the number of items in the scene.
  virtual int numberOfEntries() const = 0;
  //!
  //! @returns the item with the specified index.
  virtual CGAL::Three::Scene_item* item(Item_id) const = 0;
  //!
  //! @returns the id of the specified item.
  virtual Item_id item_id(CGAL::Three::Scene_item*) const = 0;
  //!
  //!@returns the currently selected item's index. If several items are selected, returns the last one's index.
  virtual Item_id mainSelectionIndex() const = 0;
  //!
  //!@returns the list of currently selected items indices.
  virtual QList<Item_id> selectionIndices() const = 0;
  //!
  //!@returns the index of the Item_A
  virtual Item_id selectionAindex() const = 0;
  //!
  //!@returns the index of the Item_B
  virtual Item_id selectionBindex() const = 0;

  //!
  //!@returns the scene bounding box
  virtual Bbox bbox() const = 0;
  //!
  //!@returns the length of the bounding box's diagonal.
  virtual double len_diagonal() const = 0;

public:
  /*! Emits a dataChanged signal. */
  virtual void itemChanged(Item_id i) = 0; 
  /*! Emits a dataChanged signal.*/
  virtual void itemChanged(CGAL::Three::Scene_item*) = 0;

  /*! Selects the item with the specified index. */
  virtual void setSelectedItem(Item_id) = 0;
  /*!  Adds group to the Scene. If items were selected at the
   * moment this function is called, they are added to the group.
   * If all the selected items were in the same group, the new group
   * is added as a sub-group of this group.*/
  virtual void addGroup(Scene_group_item* group) = 0;
  
}; // end interface Scene_interface
}
}

#endif // SCENE_INTERFACE_H
