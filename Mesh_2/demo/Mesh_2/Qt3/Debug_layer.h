// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_DEBUG_LAYER_H
#define CGAL_DEBUG_LAYER_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <iostream>

namespace CGAL {

class Debug_layer : public Qt_widget_layer
{
  std::ostream& stream;

public:

  Debug_layer(std::ostream& s = std::cerr) : stream(s)
  {
  }

  void draw()
  {
    stream << "redraw()" << std::endl;
  }
}; // end class Debug_layer

} // end namespace CGAL

#endif // CGAL_DEBUG_LAYER_H
