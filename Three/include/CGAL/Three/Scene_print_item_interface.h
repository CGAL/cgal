// Copyright (c) 2009,2010,2012,2015  GeometryFactory Sarl (France)
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
  virtual bool printVertexIds(CGAL::Three::Viewer_interface*) const= 0;
  //! Prints all the edges ids if their number is not too high. The limit is
  //! editable in the View menu of the application.
  //! \returns `false` if the number of ids is too high to be displayed.
  virtual bool printEdgeIds(CGAL::Three::Viewer_interface*) const= 0;
  //! Prints all the faces ids if their number is not too high. The limit is
  //! editable in the View menu of the application.
  //! \returns `false` if the number of ids is too high to be displayed.
  virtual bool printFaceIds(CGAL::Three::Viewer_interface*) const= 0;
  //! Prints all the primitive ids if their number is not too high. The limit is
  //! editable in the View menu of the application.
  virtual void printAllIds(CGAL::Three::Viewer_interface*) = 0;
  //! \brief Tests if an id should be displayed or not.
  //!
  //! \returns true if the Id should be displayed
  //! \returns false if the Id should not be displayed (if it is hidden for example)
  virtual bool testDisplayId(double, double, double, CGAL::Three::Viewer_interface*)const = 0;

  //! \brief Tests if this item should display its ids.
  //!
  //! The default behavior is to only display ids of the currently selected item (\see mainSelectionIndex()).
  //! This function allows to override this behavior.
  //! @param tets_item the currently tested TextListItem.
  //! \return true if this item should display its ids when `test_item` is tested.
  virtual bool shouldDisplayIds(CGAL::Three::Scene_item* test_item)const = 0;
};
}
}

Q_DECLARE_INTERFACE(CGAL::Three::Scene_print_item_interface, "com.geometryfactory.PolyhedronDemo.PrintInterface/1.0")
#endif // SCENE_PRINT_ITEM_INTERFACE_H

