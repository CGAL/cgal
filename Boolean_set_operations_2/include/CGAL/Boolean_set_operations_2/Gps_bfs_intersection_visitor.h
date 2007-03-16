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
// $Id$
// 
//

// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
#ifndef CGAL_GPS_BFS_INTERSECTION_VISITOR_H
#define CGAL_GPS_BFS_INTERSECTION_VISITOR_H

#include <CGAL/Boolean_set_operations_2/Gps_bfs_base_visitor.h>

CGAL_BEGIN_NAMESPACE

template <class Arrangement_>
class Gps_bfs_intersection_visitor : public Gps_bfs_base_visitor<Arrangement_>
{
  typedef  Arrangement_                                  Arrangement;
  typedef typename Arrangement::Face_iterator            Face_iterator;
  typedef typename Arrangement::Halfedge_iterator        Halfedge_iterator;
  typedef Gps_bfs_base_visitor<Arrangement>              Base;
  typedef typename Base::Edges_hash                      Edges_hash;
  typedef typename Base::Faces_hash                      Faces_hash;
 

public:

  Gps_bfs_intersection_visitor(Edges_hash* edges_hash,
                               Faces_hash* faces_hash,
                               unsigned int n_polygons): 
    Base(edges_hash, faces_hash, n_polygons)
  {}


  void flip_face(Face_iterator f1, Face_iterator f2, Halfedge_iterator he)
  {
    unsigned int ic_f2;
    ic_f2 = this->compute_ic(f1, f2, he);
    (*(this->m_faces_hash))[f2] = ic_f2;
      
    CGAL_assertion(ic_f2 <= this->m_num_of_polygons);

    // only faces that have inside counter equal to the number of polygons
    // which are intersectd, will be marked true (containted)
    if(ic_f2 == this->m_num_of_polygons)
      f2->set_contained(true);
  }

  // mark the unbounded_face (true iff contained)
  void visit_ubf(Face_iterator ubf, unsigned int ubf_ic)
  {
    CGAL_assertion(ubf->is_unbounded());
    if(ubf_ic == this->m_num_of_polygons)
      ubf->set_contained(true);
  }

  void after_scan(Arrangement&)
  {}


};
CGAL_END_NAMESPACE

#endif
