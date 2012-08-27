// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_QT_WIDGET_HYPERBOLA_2_H
#define CGAL_QT_WIDGET_HYPERBOLA_2_H

#include <CGAL/Hyperbola_2.h>
#include <CGAL/Hyperbola_ray_2.h>
#include <CGAL/IO/Qt_widget.h>

namespace CGAL {

template< class Gt >
inline
Qt_widget& operator<<(Qt_widget& s, const Hyperbola_2< Gt > &H)
{
  H.draw_qt(s);
  return s;
}

template< class Gt >
inline
Qt_widget& operator<<(Qt_widget &s, Hyperbola_ray_2<Gt> &H)
{
  H.draw_qt(s);
  return s;
}

} //namespace CGAL

#endif // CGAL_QT_WIDGET_HYPERBOLA_2_H
