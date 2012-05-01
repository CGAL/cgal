// Copyright (c) 2004,2005  INRIA Sophia-Antipolis (France).
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

#ifndef SDG_INSERT_H
#define SDG_INSERT_H


template<class SDG, class Point>
typename SDG::Vertex_handle
insert_point(SDG& sdg, const Point& p)
{
  return sdg.insert(p);
}

template<class SDG, class Point>
typename SDG::Vertex_handle
insert_segment(SDG& sdg, const Point& p1, const Point& p2)
{
  return sdg.insert(p1, p2);
}

template<class SDG, class Point>
typename SDG::Vertex_handle
insert_segment(SDG& sdg, const Point& p1, const Point& p2,
	       typename SDG::Vertex_handle v)
{
  return sdg.insert(p1, p2, v);
}


template<class SDG, class Polygon>
typename SDG::Vertex_handle
insert_polygon(SDG& sdg, const Polygon& pgn)
{
  typename SDG::Vertex_handle v;
  int psize = pgn.size();
  for (int i = 0; i < psize; i++ ) {
    v = sdg.insert( pgn[i], pgn[(i+1)%psize] );
  }
  return v;
}



#endif // SDG_INSERT_H
