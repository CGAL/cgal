// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
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
// $Id$ $Date$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef GPS_BPS_BASE_VISITOR_H
#define GPS_BPS_BASE_VISITOR_H

#include <CGAL/Unique_hash_map.h> 

CGAL_BEGIN_NAMESPACE

template <class Arrangement_>
class Gps_bfs_base_visitor
{
  typedef  Arrangement_                                  Arrangement;
  typedef typename Arrangement::Face_iterator            Face_iterator;
  typedef typename Arrangement::Halfedge_iterator        Halfedge_iterator;
public:
  typedef Unique_hash_map<Halfedge_iterator, unsigned int> Edges_hash;
  typedef Unique_hash_map<Face_iterator, unsigned int>     Faces_hash;

protected:
  Edges_hash*    m_edges_hash;
  Faces_hash*    m_faces_hash;
  unsigned int   m_num_of_polygons; // number of polygons

public:

  Gps_bfs_base_visitor(Edges_hash* edges_hash,
                       Faces_hash* faces_hash,
                       unsigned int n_pgn): 
    m_edges_hash(edges_hash),
    m_faces_hash(faces_hash),
    m_num_of_polygons(n_pgn)
  {}


  void flip_face(Face_iterator f1, Face_iterator f2, Halfedge_iterator he)
  {
    CGAL_assertion(m_edges_hash->is_defined(he) && 
                   m_edges_hash->is_defined(he->twin()) &&
                   m_faces_hash->is_defined(f1) &&
                   !m_faces_hash->is_defined(f2));

    // IC of f2 (inside counter)
    unsigned int ic_f2 = 
      (*m_faces_hash)[f1] - (*m_edges_hash)[he] + (*m_edges_hash)[he->twin()];
    (*m_faces_hash)[f2] = ic_f2;
  }

protected:

  // compute the inside count of a face
  unsigned int compute_ic(Face_iterator f1, Face_iterator f2, Halfedge_iterator he)
  {
    CGAL_assertion(m_edges_hash->is_defined(he) && 
                   m_edges_hash->is_defined(he->twin()) &&
                   m_faces_hash->is_defined(f1) &&
                   !m_faces_hash->is_defined(f2));
    unsigned int ic_f2 = 
      (*m_faces_hash)[f1] - (*m_edges_hash)[he] + (*m_edges_hash)[he->twin()];
    (*m_faces_hash)[f2] = ic_f2;

    return (ic_f2);
  }
};

CGAL_END_NAMESPACE

#endif
