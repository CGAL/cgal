// Copyright (c) 20017  GeometryFactory Sarl (France)
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
#ifndef SCENE_FACEGRAPH_INTERFACE_H
#define SCENE_FACEGRAPH_INTERFACE_H

namespace CGAL
{
namespace Three {
//! Base class for a Scene_item containing a FaceGraph
template <typename Mesh>
class Scene_facegraph_interface_item {
public:
  virtual ~Scene_facegraph_interface_item(){}
 //!Returns the item's facegraph
 virtual Mesh* facegraph() = 0;
  virtual const Mesh* facegraph()const = 0;
};
}
}
#endif // SCENE_FACEGRAPH_INTERFACE_H
