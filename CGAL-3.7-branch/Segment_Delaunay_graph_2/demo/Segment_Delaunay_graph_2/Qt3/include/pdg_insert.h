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
// $URL$
// $Id$
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>


#ifndef PDG_INSERT_H
#define PDG_INSERT_H


int get_new_id()
{
  static int gen = 0;
  gen++;

  CGAL_assertion( gen != 0 ); // guards against unwanted phenomena due to
                      // overflow
  return gen;
}



template<class PDG, class Point>
typename PDG::Vertex_handle
insert_point(PDG& pdg, const Point& p)
{
  typename PDG::Vertex_handle v = pdg.insert(p);
  v->set_info( get_new_id() );
  return v;
}

template<class PDG, class Point>
typename PDG::Vertex_handle
insert_segment(PDG& pdg, const Point& p1, const Point_2& p2, int id)
{
  typedef typename PDG::Vertex_handle Vertex_handle;

  Vertex_handle v1 = pdg.insert(p1);
  Vertex_handle v2 = pdg.insert(p2);
  Vertex_handle v3 = pdg.insert(p1,p2);

  if ( v3 == Vertex_handle() ) {
    return v3;
  }

  v1->set_info( id );
  v2->set_info( id );
  v3->set_info( id );
  return v3;
}

template<class PDG, class Point>
typename PDG::Vertex_handle
insert_segment(PDG& pdg, const Point& p1, const Point_2& p2,
	       typename PDG::Vertex_handle v, int id)
{
  typedef typename PDG::Vertex_handle Vertex_handle;

  Vertex_handle v1 = pdg.insert(p1, v);
  Vertex_handle v2 = pdg.insert(p2, v2);
  Vertex_handle v3 = pdg.insert(p1, p2, v1);

  if ( v3 == Vertex_handle() ) {
    return v3;
  }

  v1->set_info( id );
  v2->set_info( id );
  v3->set_info( id );
  return v3;
}


template<class PDG, class Point>
typename PDG::Vertex_handle
insert_segment(PDG& pdg, const Point& p1, const Point& p2)
{
  int id = get_new_id();
  return insert_segment(pdg, p1, p2, id);
}

template<class PDG, class Point>
typename PDG::Vertex_handle
insert_segment(PDG& pdg, const Point& p1, const Point& p2,
	       typename PDG::Vertex_handle v)
{
  int id = get_new_id();
  return insert_segment(pdg, p1, p2, v, id);
}


template<class PDG, class Polygon>
typename PDG::Vertex_handle
insert_polygon(PDG& pdg, const Polygon& pgn)
{
  typedef typename PDG::Vertex_handle Vertex_handle;

  int id = get_new_id();

  Vertex_handle v;
  int psize = pgn.size();
  for (int i = 0; i < psize; i++ ) {
    v = insert_segment( pdg, pgn[i], pgn[(i+1)%psize], id );
    if ( v == Vertex_handle() ) { break; }
  }
  return v;
}



#endif // PDG_INSERT_H
