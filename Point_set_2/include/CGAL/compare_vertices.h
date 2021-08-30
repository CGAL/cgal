// Copyright (c) 1999
// Max-Planck-Institute Saarbruecken (Germany). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Matthias Baesken

#ifndef CGAL_COMPARE_VERTICES_H
#define CGAL_COMPARE_VERTICES_H

#include <CGAL/license/Point_set_2.h>


namespace CGAL {

namespace internal {

// compare function objects for the priority queues used in nearest neighbor search
template<class VP, class NT,class MAP_TYPE>
class compare_vertices {
 public:
  //std::map<VP,NT,std::less<VP> > *pmap;
  MAP_TYPE* pmap;

  compare_vertices(MAP_TYPE *p){ pmap=p; }

  bool operator()(VP e1, VP e2)
  // get the priorities from the map and return result of comparison ...
  { NT& v1 = (*pmap)[e1];
    NT& v2 = (*pmap)[e2];
    return (v1 > v2);
  }
};


} // namespace internal


} //namespace CGAL

#endif // CGAL_COMPARE_VERTICES_H
