// Copyright (c) 2004,2005  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>

#ifndef SVD_INSERT_H
#define SVD_INSERT_H


template<class SVD, class Point>
typename SVD::Vertex_handle
insert_point(SVD& svd, const Point& p)
{
  return svd.insert(p);
}

template<class SVD, class Point>
typename SVD::Vertex_handle
insert_segment(SVD& svd, const Point& p1, const Point& p2)
{
  return svd.insert(p1, p2);
}

template<class SVD, class Point>
typename SVD::Vertex_handle
insert_segment(SVD& svd, const Point& p1, const Point& p2,
	       typename SVD::Vertex_handle v)
{
  return svd.insert(p1, p2, v);
}


template<class SVD, class Polygon>
typename SVD::Vertex_handle
insert_polygon(SVD& svd, const Polygon& pgn)
{
  typename SVD::Vertex_handle v;
  int psize = pgn.size();
  for (int i = 0; i < psize; i++ ) {
    v = svd.insert( pgn[i], pgn[(i+1)%psize] );
  }
  return v;
}



#endif // SVD_INSERT_H
