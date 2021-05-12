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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ophir Setter    <ophir.setter@cs.tau.ac.il>

#ifndef CGAL_GPS_BFS_JOIN_VISITOR_H
#define CGAL_GPS_BFS_JOIN_VISITOR_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/Boolean_set_operations_2/Gps_bfs_base_visitor.h>

namespace CGAL {

template <class Arrangement_>
class Gps_bfs_join_visitor :
public Gps_bfs_base_visitor<Arrangement_, Gps_bfs_join_visitor<Arrangement_> >
{
  typedef  Arrangement_                                  Arrangement;
  typedef typename Arrangement::Face_iterator            Face_iterator;
  typedef typename Arrangement::Halfedge_iterator        Halfedge_iterator;
  typedef Gps_bfs_join_visitor<Arrangement>              Self;
  typedef Gps_bfs_base_visitor<Arrangement, Self>        Base;
  typedef typename Base::Edges_hash                      Edges_hash;
  typedef typename Base::Faces_hash                      Faces_hash;

public:

  Gps_bfs_join_visitor(Edges_hash* edges_hash, Faces_hash* faces_hash, unsigned int n_pgn):
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
    // at least one polygon contains the face.
    return (ic > 0);
  }

  void after_scan(Arrangement&)
  {}

};

} //namespace CGAL

#endif
