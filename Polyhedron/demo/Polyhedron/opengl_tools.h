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
// 
//
// Author(s)     : Laurent Rineau <laurent.rineau__CGAL@normalesup.org>

#ifndef CGAL_OPENGL_TOOLS_H
#define CGAL_OPENGL_TOOLS_H

//#ifdef CGAL_GLEW_ENABLED
//#else
# include <CGAL/gl.h>
//#endif

namespace CGAL {
namespace GL {

class Color {
  GLfloat c[4];
public:
  Color() {
    ::glGetFloatv(GL_CURRENT_COLOR, &c[0]);
  }
  ~Color() {
    set_rgb_color(c[0], c[1], c[2], c[3]);
  }
  void set_rgb_color(GLfloat r, GLfloat g, GLfloat b, GLfloat a = 1.f) {
    ::glColor4f(r, g, b, a);
  }
}; // end class Color;

class Point_size {
  GLfloat ps;
public:
  Point_size() {
    ::glGetFloatv(GL_POINT_SIZE, &ps);
  }
  ~Point_size() {
    set_point_size(ps);
  }
  void set_point_size(GLfloat v) {
    ::glPointSize(v);
  }
}; // end class Point_size

} // end namespace GL
} // end namespace CGAL

#endif // not CGAL_OPENGL_TOOLS_H
