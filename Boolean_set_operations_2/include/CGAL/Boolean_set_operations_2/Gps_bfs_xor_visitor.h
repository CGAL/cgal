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

#ifndef CGAL_GPS_BFS_XOR_VISITOR_H
#define CGAL_GPS_BFS_XOR_VISITOR_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/Boolean_set_operations_2/Gps_bfs_base_visitor.h>

namespace CGAL {

template <typename Arrangement_>
class Gps_bfs_xor_visitor :
    public Gps_bfs_base_visitor<Arrangement_, Gps_bfs_xor_visitor<Arrangement_>> {
  using Arrangement = Arrangement_;
  using Face_iterator = typename Arrangement::Face_iterator;
  using Halfedge_iterator = typename Arrangement::Halfedge_iterator;
  using Self = Gps_bfs_xor_visitor<Arrangement>;
  using Base = Gps_bfs_base_visitor<Arrangement, Self>;
  using Edges_hash = typename Base::Edges_hash;
  using Faces_hash = typename Base::Faces_hash;

public:
  Gps_bfs_xor_visitor(Edges_hash* edges_hash, Faces_hash* faces_hash,
                      std::size_t n_pgn) :
    Base(edges_hash, faces_hash, n_pgn)
  {}

  //! contained_criteria
/*! contained_criteria is used to the determine if the face which has
  inside count should be marked as contained.
  \param ic the inner count of the talked-about face.
  \return true if the face of ic, otherwise false.
*/
  bool contained_criteria(std::size_t ic) {
    // xor means odd number of polygons.
    return (ic % 2) == 1;
  }

  //! after_scan postprocessing after bfs scan.
  /*! The function fixes some of the curves, to be in the same direction as the
   * half-edges.
   *
   * \param arr The given arrangement.
   */
  void after_scan(Arrangement& arr) {
    using Traits = typename Arrangement::Geometry_traits_2;
    using X_monotone_curve_2 = typedef typename Traits::X_monotone_curve_2;

    Traits tr;
    auto cmp_endpoints = tr.compare_endpoints_xy_2_object();
    auto ctr_opp = tr.construct_opposite_2_object();

    for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
      Halfedge_iterator he = eit;
      const X_monotone_curve_2& cv = he->curve();
      const bool is_cont = he->face()->contained();
      const Comparison_result he_res =
        ((Arr_halfedge_direction)he->direction() == ARR_LEFT_TO_RIGHT) ?
        SMALLER : LARGER;
      const bool has_same_dir = (cmp_endpoints(cv) == he_res);

      if ((is_cont && !has_same_dir) || (!is_cont && has_same_dir))
        arr.modify_edge(he, ctr_opp(cv));
    }
  }
};

} //namespace CGAL

#endif
