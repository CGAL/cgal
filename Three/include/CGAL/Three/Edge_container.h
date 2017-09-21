// Copyright (c) 2017  GeometryFactory Sarl (France)
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
// Author(s)     : Maxime Gimeno

#ifndef EDGE_CONTAINER_H
#define EDGE_CONTAINER_H

#include <CGAL/license/Three.h>


#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Primitive_container.h>

using namespace CGAL::Three;
#ifdef demo_framework_EXPORTS
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_EXPORT
#else
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_IMPORT
#endif
namespace CGAL {
namespace Three {

struct DEMO_FRAMEWORK_EXPORT Edge_container :public Primitive_container
{

  enum vbosName {
    Vertices = 0,
    Indices,
    NbOfVbos
  };

    Edge_container(int program, bool indexed);
    void initGL(Scene_item *item, Viewer_interface *viewer) const Q_DECL_OVERRIDE;
    void draw(const Scene_item &item, CGAL::Three::Viewer_interface* viewer,
              bool is_color_uniform = true) const Q_DECL_OVERRIDE;
    //drawing variables
    QVector4D plane;
}; //end of class Triangle_container
}
}

#endif // EDGE_CONTAINER_H
