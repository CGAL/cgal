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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Laurent RINEAU

//! \file Scene_interface.h


#ifndef SCENE_INTERFACE_H
#define SCENE_INTERFACE_H
#include <CGAL/license/Three.h>
/*!
* \ingroup PkgThree
* The RenderingMode determines which of an item's primitives must be drawn.
* It can be Points, PointsPlusNormals, Wireframe, Flat, FlatPlusEdges, or Gouraud.
* - Points, PointsPlusNormals, and Wireframe have no light model.
* - Flat and FlatPlusEdges use a basic light model with one normal per facet.
* - Gouraud uses the same light model but with one normal per vertex.
*/
enum RenderingMode
{
  Points = 0, //! Renders only points without lighting.
  PointsPlusNormals, //!Renders points and normals.
  Wireframe, //!Renders only edges.
  Flat, //!Renders only faces, with a lighting per face.
  FlatPlusEdges, //!Renders flat faces and edges.
  Gouraud, //!Renders only faces, with a lighting per vertex.
  ShadedPoints, //!Renders only points with lighting.
  NumberOfRenderingMode //!Number of values in this enum.
};


#include <QString>
#include <QColor>
#include <QList>
#include <algorithm>
#include <cmath>
#include <CGAL/Bbox_3.h>
namespace CGAL{
namespace Three{
class Scene_item;
class Scene_group_item;
}
}

namespace CGAL {
namespace Three{
/*!
 * This is the class given to the plugins to interact with the scene.
 * @todo Laurent Rineau Scene_interface should be a QObject, with the useful signals and slots of the scene

 * */
class Scene_interface {
public:

  //!A bounding box is a box with each face corresponding to an extremum of its contents.

  typedef CGAL::Bbox_3 Bbox;

  //!\brief Integer used as the index of a Scene_item.
  //! An item's index is its position in the Geometric Objects list.
  typedef int Item_id;
  virtual ~Scene_interface() {};
  //! Adds `item` to the Geometric Objects list.
  //! The scene then takes ownership of `item`.
  //! If `update_bbox` is `true`, the whole bbox of the scene is re-computed,
  //! along with the global offset.
  //!@returns the index of the new item.
  virtual Item_id addItem(CGAL::Three::Scene_item* item) = 0;
  //! \brief Replaces an item by a new one in the scene.
  //! The item which id is `id` is replaced by `item`.
  //! The first one is destroyed and gives its index to the second one.
  //!@returns a pointer to the old item.
  virtual Scene_item* replaceItem(Item_id id, CGAL::Three::Scene_item* item) = 0;
  //!Moves item to the targeted group.
  virtual void changeGroup(CGAL::Three::Scene_item* item, CGAL::Three::Scene_group_item* target_group) = 0;

  /*! Erases an item in the list and destroy it.
   * @returns the index of the item just before the one that is erased,
   * or just after.
   * @returns -1 if the list is empty.*/
  virtual Item_id erase(Item_id) = 0;
  /*! Erase the items with the target indices and destroy them.
   * @returns the index of the polyhedron just before the
   * one that is erased, or just after. Returns -1 if
   * the list is empty.
   */
  virtual int erase(QList<int>) = 0;

  /*! Creates a copy of the `id`th item .
   * @returns the index of the new item (-1 on error).
   */
    virtual Item_id duplicate(Item_id id) = 0;

  // Accessors (getters)
  //! \brief The total number of items in the scene.
  //!@returns the number of items in the scene.
  virtual int numberOfEntries() const = 0;
  //!\brief The `id`th item.
  //! @returns the item with the specified index.
  virtual CGAL::Three::Scene_item* item(Item_id id) const = 0;
  //!\brief The id of `item`
  //! @returns the id of the specified item.
  virtual Item_id itemId(CGAL::Three::Scene_item* item) const = 0;
  //!\brief The currently selected item's index.
  //!@returns the currently selected item's index.
  //!@returns -1 if none or several items are selected
  virtual Item_id mainSelectionIndex() const = 0;
  //!The id of the currently selected item.
  //!@returns the list of currently selected items indices.
  virtual QList<Item_id> selectionIndices() const = 0;

  //!\brief The scene's Bbox
  //!@returns the scene's bounding box
  //! @see Scene_interface::Bbox
  virtual Bbox bbox() const = 0;
  //!\brief The Bbox of the visible items of the Bbox.
  //!@returns the scene's bounding box, computed only for currently visible items.
  //! @see Scene_interface::Bbox
  virtual Bbox visibleBbox() const = 0;
  //!The length of the diagonal of the scene's Bbox
  //!@returns the length of the bounding box's diagonal.
  virtual double bboxDiagonalLength() const = 0;

public:
  //! Updates the information about the `i`th item in the
  //! Geometric Objects list and redraws the scene.
  //!
  //! This should be called when the tooltip, graphical Tooltip,
  //! visibility, name or color of an existing item is changed.
  virtual void itemChanged(Item_id i) = 0;
  //! Updates the information about `item` in the
  //! Geometric Objects list and redraws the scene.
  //!
  //! This should be called when the tooltip, graphical Tooltip,
  //! visibility, name or color of an existing item is changed.
  virtual void itemChanged(CGAL::Three::Scene_item* item) = 0;
  //! Called when the visibility of an item has changed to re compute
  //! the scene Bbox without recentering it. This allows to keep an
  //! adapted frustum without moving the camera.
  virtual void itemVisibilityChanged(CGAL::Three::Scene_item*) = 0;
  //! Clears the current selection then sets the selected item to the target index.
  //! Used to update the selection in the Geometric Objects view.
  virtual void setSelectedItem(Item_id) = 0;
  //! Sets the target list of indices as the selected indices and returns it.
  virtual QList<int> setSelectedItemsList(QList<int> l ) = 0;
}; // end interface Scene_interface
}
}


#endif // SCENE_INTERFACE_H
