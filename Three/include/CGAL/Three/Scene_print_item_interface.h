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


//! An item that wants to print its primitive IDs must derive from this interface.
class Scene_print_item_interface {
public:
  virtual ~Scene_print_item_interface(){}
 //! Finds the spot the closest to `point` and prints the id of the corresponding Primitive (vertex, edge or face).
 virtual void printPrimitiveId(QPoint, CGAL::Three::Viewer_interface*) = 0;
  //! Prints all the primitive ids if their number is not too high. The limit is
  //! editable in the View menu of the application.
 virtual void printPrimitiveIds(CGAL::Three::Viewer_interface*)const = 0;
  //! \brief Tests if an id should be displayed or not.
  //!
  //! \returns true if the Id should be displayed
  //! \returns false if the Id should not be displayed (if it is hidden for example)
  virtual bool testDisplayId(double, double, double, CGAL::Three::Viewer_interface*)const = 0;
};
}
}

Q_DECLARE_INTERFACE(CGAL::Three::Scene_print_item_interface, "com.geometryfactory.PolyhedronDemo.PrintInterface/1.0")
#endif // SCENE_PRINT_ITEM_INTERFACE_H

