// Copyright (c) 2003  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Laurent Rineau <laurent.rineau__CGAL@normalesup.org>

#ifndef CGAL_OPENGL_TOOLS_H
#define CGAL_OPENGL_TOOLS_H

#include <CGAL/Three/Viewer_interface.h>

namespace CGAL {
namespace GL {

class Point_size {
  GLfloat ps;
  CGAL::Three::Viewer_interface* viewer;
public:
  Point_size(CGAL::Three::Viewer_interface* viewer) : viewer(viewer) {
    viewer->glGetFloatv(GL_POINT_SIZE, &ps);
  }
  ~Point_size() {
    set_point_size(ps);
  }
  void set_point_size(GLfloat v) {
    viewer->glPointSize(v);
  }
}; // end class Point_size

} // end namespace GL
} // end namespace CGAL

#endif // not CGAL_OPENGL_TOOLS_H
