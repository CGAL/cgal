// Copyright (c) 2008  GeometryFactory, Sophia-Antipolis (France).
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

#ifndef CGAL_MESH_3_DIHEDRAL_ANGLE_3_H
#define CGAL_MESH_3_DIHEDRAL_ANGLE_3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/number_type_basic.h>
#include <CGAL/Kernel_traits.h>
#include <cmath>


namespace CGAL {
namespace Mesh_3 {

/**
 * Computes dihedral angle of planes (a,b,c) and (a,b,d)
 */
template <typename K>
typename K::FT
dihedral_angle(const typename K::Point_3& a,
               const typename K::Point_3& b,
               const typename K::Point_3& c,
               const typename K::Point_3& d,
               K k = K())
{
  // Now in the CGAL kernels
  return k.compute_approximate_dihedral_angle_3_object()(a, b, c, d);
}


/**
 * Computes dihedral angle of planes (a,b,c) and (a,b,d)
 */
template <typename Point_3>
typename Kernel_traits<Point_3>::Kernel::FT
dihedral_angle(const Point_3& a, const Point_3& b,
               const Point_3& c, const Point_3& d)
{
  return
    CGAL::Mesh_3::dihedral_angle(a, b, c, d,
                                 typename Kernel_traits<Point_3>::Kernel());
}

} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_DIHEDRAL_ANGLE_3_H
