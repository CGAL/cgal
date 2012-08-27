// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// $Date$
// 
//
// Author(s)     : Ron Wein  <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_QT_WIDGET_LINEAR_OBJECT_2_H
#define CGAL_QT_WIDGET_LINEAR_OBJECT_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Arr_linear_traits_2.h>

namespace CGAL {

/*!
 * Export a polyline to a window stream 
 */
template <class K>
Qt_widget & operator<<(Qt_widget & ws, const Arr_linear_object_2<K> & o)
{
  if(o.is_segment())
  {
    ws << o.segment();
    return ws;
  }
  if(o.is_ray())
  {
    ws << o.ray();
    return ws;
  }

  CGAL_assertion(o.is_line());
  ws << o.line();
  return ws;
}

} //namespace CGAL

#endif
