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


#ifndef CGAL_SURFACE_MESHER_NULL_ORACLE_VISITOR_H
#define CGAL_SURFACE_MESHER_NULL_ORACLE_VISITOR_H

#include <CGAL/license/Surface_mesher.h>


namespace CGAL {

  namespace Surface_mesher {

  /** \interface OracleVisitor
      \brief Concept of a visitor of oracles.

      \fn void new_point(Point& p)
      Called before a new Point is returned by the oracle.
      \param p The point that has just been computed and will be returned.
  */

  /** Trivial model of the OracleVisitor concept. */
  struct Null_oracle_visitor
  {
    template <class P>
    void new_point(P&) const
    {
    }
  };

  }  // namespace Surface_mesher

} // namespace CGAL


#endif  // CGAL_SURFACE_MESHER_NULL_ORACLE_VISITOR_H
