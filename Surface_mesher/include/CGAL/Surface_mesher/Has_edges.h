// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_SURFACE_MESHER_HAS_EDGES_H
#define CGAL_SURFACE_MESHER_HAS_EDGES_H

#include <CGAL/license/Surface_mesher.h>


namespace CGAL {
  namespace Surface_mesher {

    struct Has_edges {
      static bool has_edges() { return true; }
    };
    struct Has_no_edges {
      static bool has_edges() { return false; }
    };

  } // end namespace Surface_mesher
} // end namespace CGAL

#endif // CGAL_SURFACE_MESHER_HAS_EDGES_H
