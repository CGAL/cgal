// Copyright (c) 2008 ETH Zurich (Switzerland)
// Copyright (c) 2008-2009  INRIA Sophia-Antipolis (France)
// Copyright (c) 2009  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Laurent Rineau, Camille Wormser, Jane Tournois, Pierre Alliez

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_BBOX_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_BBOX_3_DO_INTERSECT_H

// Turn off Visual C++ warning
#ifdef _MSC_VER
#pragma warning ( disable : 4003 )
#endif

#include <CGAL/Bbox_3.h>

namespace CGAL {
  bool
  inline
  do_intersect(const CGAL::Bbox_3& c,
               const CGAL::Bbox_3& bbox)
  {
    return CGAL::do_overlap(c, bbox);
  }

} //namespace CGAL

#endif  // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_BBOX_3_DO_INTERSECT_H
