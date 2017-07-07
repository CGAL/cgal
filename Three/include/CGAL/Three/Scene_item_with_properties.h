// Copyright (c) 2016  GeometryFactory Sarl (France)
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

#ifndef SCENE_ITEM_WITH_PROPERTIES_H
#define SCENE_ITEM_WITH_PROPERTIES_H

#include <CGAL/license/Three.h>

namespace CGAL
{
namespace Three {
  class Scene_item;

//! Base class to allow an item to copy properties from another.
//! Properties reprensent the current state of an item : its color,
//! the position of its manipulated frame, ...
class Scene_item_with_properties {
public:
  virtual ~Scene_item_with_properties(){}
 //!\brief Copies properties from another Scene_item.
 //!
 //! Override this function to specify what must be copied.
 virtual void copyProperties(Scene_item*)=0; //pure virtual method
};
}
}
#endif // SCENE_ITEM_WITH_PROPERTIES_H

