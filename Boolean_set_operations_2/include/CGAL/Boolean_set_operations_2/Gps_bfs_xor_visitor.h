// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ophir Setter    <ophir.setter@cs.tau.ac.il>

#ifndef CGAL_GPS_BFS_XOR_VISITOR_H
#define CGAL_GPS_BFS_XOR_VISITOR_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/Boolean_set_operations_2/Gps_bfs_base_visitor.h>

namespace CGAL {

template <class Arrangement_>
class Gps_bfs_xor_visitor : 
public Gps_bfs_base_visitor<Arrangement_, Gps_bfs_xor_visitor<Arrangement_> >
{
  typedef  Arrangement_                                  Arrangement;
  typedef typename Arrangement::Face_iterator            Face_iterator;
  typedef typename Arrangement::Halfedge_iterator        Halfedge_iterator;
  typedef Gps_bfs_xor_visitor<Arrangement>               Self;
  typedef Gps_bfs_base_visitor<Arrangement, Self>        Base;
  typedef typename Base::Edges_hash                      Edges_hash;
  typedef typename Base::Faces_hash                      Faces_hash;

public:

  Gps_bfs_xor_visitor(Edges_hash* edges_hash, Faces_hash* faces_hash, 
                      unsigned int n_pgn) :
    Base(edges_hash, faces_hash, n_pgn)
  {}

    //! contained_criteria
/*! contained_criteria is used to the determine if the face which has 
  inside count should be marked as contained.
  \param ic the inner count of the talked-about face.
  \return true if the face of ic, otherwise false.
*/
  bool contained_criteria(unsigned int ic)
  {
    // xor means odd number of polygons.
    return (ic % 2) == 1;
  }

  //! after_scan post-processing after bfs scan.
/*! The function fixes some of the curves, to be in the same direction as the
    half-edges. 
  
  \param arr The given arrangment.
*/
  void after_scan(Arrangement& arr)
  {
    typedef typename Arrangement::Geometry_traits_2 Traits;
    typedef typename Traits::Compare_endpoints_xy_2 Compare_endpoints_xy_2;
    typedef typename Traits::Construct_opposite_2   Construct_opposite_2;
    typedef typename Traits::X_monotone_curve_2     X_monotone_curve_2;
    typedef typename Arrangement::Edge_iterator     Edge_iterator;

    Traits tr;
    Compare_endpoints_xy_2 cmp_endpoints =
      tr.compare_endpoints_xy_2_object();
    Construct_opposite_2 ctr_opp = tr.construct_opposite_2_object();

    for(Edge_iterator eit = arr.edges_begin();
        eit != arr.edges_end();
        ++eit)
    {
      Halfedge_iterator         he = eit;
      const X_monotone_curve_2& cv = he->curve();
      const bool                is_cont = he->face()->contained();
      const Comparison_result   he_res = 
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
