// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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
//
//
// Author(s)     : Laurent Rineau


#ifndef CGAL_SURFACE_MESHER_POINT_SURFACE_INDICES_VISITOR_H
#define CGAL_SURFACE_MESHER_POINT_SURFACE_INDICES_VISITOR_H

#include <CGAL/license/Surface_mesher.h>


namespace CGAL {

  namespace Surface_mesher {

  /** Model of the OracleVisitor concept.
      This model of OracleVisitor sets the point "surface_index" to a
      constant \c int.
   */
  struct Point_surface_indices_visitor
  {
    int i;

    Point_surface_indices_visitor(const int index) : i(index)
    {
    }

    template <class P>
    void new_point(P& p) const
    {
      p.set_surface_index(i);
    }
  }; // end class Point_surface_indices_visitor

  }  // namespace Surface_mesher

} // namespace CGAL


#endif  // CGAL_SURFACE_MESHER_POINT_SURFACE_INDICES_VISITOR_H
