// Copyright (c) 2009,2010,2012,2015  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime GIMENO

#ifndef SCENE_PRINT_ITEM_INTERFACE_H
#define SCENE_PRINT_ITEM_INTERFACE_H

#include <CGAL/license/Three.h>
#include <QtPlugin>
#include <QPoint>
namespace CGAL
{
namespace Three {
  class Viewer_interface;
  class Scene_item;


//! An item that wants to print its primitive IDs must derive from this interface.
class Scene_print_item_interface {
public:
  virtual ~Scene_print_item_interface(){}
 //! Finds the spot the closest to `point` and prints the id of the corresponding Primitive (vertex, edge or face).
 virtual void printPrimitiveId(QPoint, CGAL::Three::Viewer_interface*) = 0;

  //! Prints all the vertices ids if their number is not too high. The limit is
  //! editable in the View menu of the application.
  //! \returns `false` if the number of ids is too high to be displayed.
  virtual bool printVertexIds() const= 0;
  //! Prints all the edges ids if their number is not too high. The limit is
  //! editable in the View menu of the application.
  //! \returns `false` if the number of ids is too high to be displayed.
  virtual bool printEdgeIds() const= 0;
  //! Prints all the faces ids if their number is not too high. The limit is
  //! editable in the View menu of the application.
  //! \returns `false` if the number of ids is too high to be displayed.
  virtual bool printFaceIds() const= 0;
  //! Prints all the primitive ids if their number is not too high. The limit is
  //! editable in the View menu of the application.
  virtual void printAllIds() = 0;
  //! \brief tests if an id should be displayed or not.
  //!
  //! \returns true if the Id should be displayed
  //! \returns false if the Id should not be displayed (if it is hidden for example)
  virtual bool testDisplayId(double, double, double, CGAL::Three::Viewer_interface*)const = 0;

  //! \brief tests if this item should display its ids.
  //!
  //! The default behavior is to only display ids of the currently selected item (\see mainSelectionIndex()).
  //! This function allows to override this behavior.
  //! @param test_item the currently tested TextListItem.
  //! \return true if this item should display its ids when `test_item` is tested.
  virtual bool shouldDisplayIds(CGAL::Three::Scene_item* test_item)const = 0;
};
}
}

Q_DECLARE_INTERFACE(CGAL::Three::Scene_print_item_interface, "com.geometryfactory.PolyhedronDemo.PrintInterface/1.0")
#endif // SCENE_PRINT_ITEM_INTERFACE_H

