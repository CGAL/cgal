// Copyright (c) 2012-2015  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

//! \file Scene_interface.h


#ifndef SCENE_INTERFACE_H
#define SCENE_INTERFACE_H
#include <CGAL/license/Three.h>
/*!
* \ingroup PkgThreeRef
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
  GouraudPlusEdges, //!Renders faces with a lighting per vertex, and edges.
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
 * */
class Scene_interface {
public:

  //!A bounding box is a box with each face corresponding to an extremum of its contents.

  typedef CGAL::Bbox_3 Bbox;

  //!\brief Integer used as the index of a Scene_item.
  //! An item's index is its position in the Geometric Objects list.
  typedef int Item_id;
  virtual ~Scene_interface() {};
  //!Adds an item to the Geometric Objects list.
  //!@returns the index of the new item.
  virtual Item_id addItem(CGAL::Three::Scene_item* item) = 0;
  //!Adds a CGAL::Three::Scene_item* to the list of children.
  virtual void addChild(Scene_item* item)=0;
  //! \brief Replaces an item by a new one in the scene.
  //! The item which id is `id` is replaced by `item`.
  //! The first one is deleted and gives its index to the second one.
  //! If emit_item_about_to_be_destroyed is true, emits
  //! an itemAboutToBeDestroyed signal.
  //!@returns a pointer to the old item.
  virtual Scene_item* replaceItem(Item_id id, CGAL::Three::Scene_item* item, bool emit_item_about_to_be_destroyed = false) = 0;
  //!Moves item to the targeted group.
  virtual void changeGroup(CGAL::Three::Scene_item* item, CGAL::Three::Scene_group_item* target_group) = 0;

  /*! Erases an item in the list.
   * @returns the index of the item just before the one that is erased,
   * or just after.
   * @returns -1 if the list is empty.*/
  virtual Item_id erase(Item_id) = 0;
  /*! Deletes the items with the target indices.
   * @returns the index of the polyhedron just before the
   * one that is erased, or just after. Returns -1 if
   * the list is empty.
   */
  virtual int erase(QList<int>) = 0;

  /*! Creates a copy of the item whith the id `id`.
   * @returns the index of the new item (-1 on error).
   */
    virtual Item_id duplicate(Item_id id) = 0;

  // Accessors (getters)
  //! \brief The number of items
  //!@returns the number of items in the scene.
  virtual int numberOfEntries() const = 0;
  //!\brief The `id`th item.
  //! @returns the item with the specified index.
  virtual CGAL::Three::Scene_item* item(Item_id id) const = 0;
  //!\brief The id of `item`
  //! @returns the id of the specified item.
  virtual Item_id item_id(CGAL::Three::Scene_item* item) const = 0;
  //!\brief The currently selected item's index.
  //!@returns the currently selected item's index.
  //!@returns -1 if none or several items are selected
  virtual Item_id mainSelectionIndex() const = 0;
  //!The id of the currently selected item.
  //!@returns the list of currently selected items indices.
  virtual QList<Item_id> selectionIndices() const = 0;
  //!Item_A is designated with the column A/B in the Geometric Objetcts widget.
  //!@returns the index of the Item_A
  virtual Item_id selectionAindex() const = 0;
  //!Item_B is designated with the column A/B in the Geometric Objetcts widget.
  //!@returns the index of the Item_B
  virtual Item_id selectionBindex() const = 0;

  //!\brief The scene's Bbox
  //!@returns the scene's bounding box
  //! @see Scene_interface::Bbox
  virtual Bbox bbox() const = 0;
  //!The length of the diagonal of the scene's Bbox
  //!@returns the length of the bounding box's diagonal.
  virtual double len_diagonal() const = 0;

public:
  //! Updates the information about the `i`th item in the
  //! Geometric Objects list and redraws the scene.
  virtual void itemChanged(Item_id i) = 0;
  //! Updates the information about `item` in the
  //! Geometric Objects list and redraws the scene.
  virtual void itemChanged(CGAL::Three::Scene_item* item) = 0;
  //! Re computes the scene Bbox without recentering it.
  virtual void itemVisibilityChanged(CGAL::Three::Scene_item*) = 0;
  //! Clears the current selection then sets the selected item to the target index.
  //! Used to update the selection in the Geometric Objects view.
  virtual void setSelectedItem(Item_id) = 0;
  //! \brief ignore data updating.
  //!
  //! This will ignore all the individual calls to `itemChanged()` until
  //! `setUpdatesEnabled()` is called whith `b` being `true`.
  //!
  virtual void setUpdatesEnabled(bool b) =0;
  //!
  //! \brief Updates all the items in the SceneView.
  //!
  virtual void allItemsChanged() = 0;
}; // end interface Scene_interface
}
}


#endif // SCENE_INTERFACE_H
