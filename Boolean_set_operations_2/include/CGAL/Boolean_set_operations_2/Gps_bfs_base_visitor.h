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
// Author(s) : Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ophir Setter    <ophir.setter@cs.tau.ac.il>

#ifndef CGAL_GPS_BPS_BASE_VISITOR_H
#define CGAL_GPS_BPS_BASE_VISITOR_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/Unique_hash_map.h>

namespace CGAL {

//! Gps_bfs_base_visitor
/*! This is a base class for all visitors that are responsible for merging
 * polygon sets.
 * We use DerivedVisitor for static polymorphism for using contained_criteria
 * which determines if we should mark the face as contained given the inside
 * count of the face.
 */
template <class Arrangement_, class DerivedVisitor>
class Gps_bfs_base_visitor {
  using Arrangement = Arrangement_;
  using Face_iterator = typename Arrangement::Face_iterator;
  using Halfedge_iterator = typename Arrangement::Halfedge_iterator;

public:
  using Edges_hash = Unique_hash_map<Halfedge_iterator, std::size_t>;
  using Faces_hash = Unique_hash_map<Face_iterator, std::size_t>;

protected:
  Edges_hash* m_edges_hash;
  Faces_hash* m_faces_hash;
  std::size_t m_num_of_polygons; // number of polygons

public:

  Gps_bfs_base_visitor(Edges_hash* edges_hash,
                       Faces_hash* faces_hash,
                       std::size_t n_pgn):
    m_edges_hash(edges_hash),
    m_faces_hash(faces_hash),
    m_num_of_polygons(n_pgn)
  {}

  //! discovered_face
  /*! discovered_face is called by Gps_bfs_scanner when it reveals a new face
   * during a BFS scan. In the BFS traversal we are going from old_face to
   * new_face through the half-edge he.
   * \param old_face The face that was already revealed
   * \param new_face The face that we have just now revealed
   * \param he The half-edge that is used to traverse between them.
   */
  void discovered_face(Face_iterator old_face,
                       Face_iterator new_face,
                       Halfedge_iterator he) {
    std::size_t ic = compute_ic(old_face, new_face, he);

    if (static_cast<DerivedVisitor*>(this)->contained_criteria(ic))
      new_face->set_contained(true);
  }

  // mark the unbounded_face (true iff contained)
  void visit_ubf(Face_iterator ubf, std::size_t ubf_ic) {
    CGAL_assertion(ubf->is_unbounded());
    if (static_cast<DerivedVisitor*>(this)->contained_criteria(ubf_ic))
      ubf->set_contained(true);
  }

protected:
  // compute the inside count of a face
  std::size_t compute_ic(Face_iterator f1,
                          Face_iterator f2,
                          Halfedge_iterator he) {
    CGAL_assertion(m_edges_hash->is_defined(he) &&
                   m_edges_hash->is_defined(he->twin()) &&
                   m_faces_hash->is_defined(f1) &&
                   !m_faces_hash->is_defined(f2));
    std::size_t ic_f2 =
      (*m_faces_hash)[f1] - (*m_edges_hash)[he] + (*m_edges_hash)[he->twin()];
    (*m_faces_hash)[f2] = ic_f2;

    return (ic_f2);
  }
};

} //namespace CGAL

#endif
