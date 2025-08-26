// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//

// Author(s):      Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ophir Setter    <ophir.setter@cs.tau.ac.il>

#ifndef CGAL_GPS_BFS_INTERSECTION_VISITOR_H
#define CGAL_GPS_BFS_INTERSECTION_VISITOR_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/Boolean_set_operations_2/Gps_bfs_base_visitor.h>

namespace CGAL {

template <class Arrangement_>
class Gps_bfs_intersection_visitor :
    public Gps_bfs_base_visitor<Arrangement_, Gps_bfs_intersection_visitor<Arrangement_>> {
  using Arrangement = Arrangement_;
  using Face_iterator = typename Arrangement::Face_iterator;
  using Halfedge_iterator = typename Arrangement::Halfedge_iterator;
  using Self = Gps_bfs_intersection_visitor<Arrangement>;
  using Base = Gps_bfs_base_visitor<Arrangement, Self>;
  using Edges_hash = typename Base::Edges_hash;
  using Faces_hash = typename Base::Faces_hash;

public:
  Gps_bfs_intersection_visitor(Edges_hash* edges_hash,
                               Faces_hash* faces_hash,
                               std::size_t n_polygons):
    Base(edges_hash, faces_hash, n_polygons)
  {}

  //! contained_criteria
  /*! contained_criteria is used to the determine if the face which has
   * inside count should be marked as contained.
   * \param ic the inner count of the talked-about face.
   * \return true if the face of ic, otherwise false.
   */
  bool contained_criteria(std::size_t ic) {
    // intersection means that all polygons contain the face.
    CGAL_assertion(ic <= this->m_num_of_polygons);
    return (ic == this->m_num_of_polygons);
  }

  void after_scan(Arrangement&) {}
};

} //namespace CGAL

#endif
