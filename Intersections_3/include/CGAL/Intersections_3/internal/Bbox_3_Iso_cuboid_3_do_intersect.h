// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// Copyright (c) 2010, 2014 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_ISO_CUBOID_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_ISO_CUBOID_3_DO_INTERSECT_H

#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {
  
namespace Intersections {

namespace internal {

  template <class K>
  bool do_intersect(const CGAL::Bbox_3& bb,
                    const typename K::Iso_cuboid_3& ic,
                    const K& k)
  {
      if (bb.xmax() < ic.xmin() || ic.xmax() < bb.xmin())
        return false;
    if (bb.ymax() < ic.ymin() || ic.ymax() < bb.ymin())
        return false;
    if (bb.zmax() < ic.zmin() || ic.zmax() < bb.zmin())
        return false;
    return true;
  }

} // namespace internal
} // namespace Intersections
} //namespace CGAL

#endif  // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_ISO_CUBOID_3_DO_INTERSECT_H
