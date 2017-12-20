// Copyright (c) 1999  
// Max-Planck-Institute Saarbruecken (Germany). All rights reserved.
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
