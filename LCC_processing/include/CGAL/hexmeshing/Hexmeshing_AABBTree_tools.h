// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>, Théo Bénard <benard320@gmail.com>
//

#ifndef CGAL_HEXMESHING_AABBTREE_TOOLS_H
#define CGAL_HEXMESHING_AABBTREE_TOOLS_H

#include <CGAL/license/LCC_processing.h>

#include <CGAL/Kernel_traits.h>
#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/Side_of_triangle_mesh.h>

namespace CGAL::internal::Hexmeshing
{
  template<typename Tree>
  bool is_intersect(double x1, double y1, double z1,
                    double x2, double y2, double z2,
                    double x3, double y3, double z3,
                    double x4, double y4, double z4,
                    const Tree& t)
  {
    using Kernel=typename CGAL::Kernel_traits<typename Tree::Point>::Kernel;
    typename Tree::Point p1(x1,y1,z1);
    typename Tree::Point p2(x2,y2,z2);
    typename Tree::Point p3(x3,y3,z3);
    typename Tree::Point p4(x4,y4,z4);

    // And compute the two triangles
    typename Kernel::Triangle_3 t1(p1, p2, p3);
    if(t.do_intersect(t1))
    { return true; }

    t1=typename Kernel::Triangle_3(p1, p3, p4);
    if(t.do_intersect(t1))
    { return true; }

    return false;
  }

  template<typename Tree>
  bool is_intersect(double x1, double y1, double z1,
                    double x2, double y2, double z2,
                    const Tree& t)
  {
    return
        is_intersect(x1,y1,z1, x2,y1,z1, x2,y1,z2, x1,y1,z2, t) || // f1 y1
        is_intersect(x2,y2,z1, x2,y1,z1, x2,y1,z2, x2,y2,z2, t) || // f2 x2
        is_intersect(x1,y2,z1, x2,y2,z1, x2,y2,z2, x1,y2,z2, t) || // f3 y2
        is_intersect(x1,y1,z1, x1,y1,z2, x1,y2,z2, x1,y2,z1, t) || // f4 x1
        is_intersect(x1,y1,z1, x2,y1,z1, x2,y2,z1, x1,y2,z1, t) || // f5 z1
        is_intersect(x1,y1,z2, x2,y1,z2, x2,y2,z2, x1,y2,z2, t);   // f6 z2
  }

  template<typename LCC, typename Tree>
  bool is_intersect(const LCC& lcc, typename LCC::Dart_const_descriptor dh, const Tree& t)
  {
    CGAL::Bbox_3 bbox=lcc.point(dh).bbox();
    // For each vertex of the volume
    for(auto it=lcc.template one_dart_per_incident_cell<0,3>(dh).begin(),
        itend=lcc.template one_dart_per_incident_cell<0,3>(dh).end(); it!=itend; ++it)
    { bbox+=lcc.point(it).bbox(); }

    return is_intersect(bbox.xmin(), bbox.ymin(), bbox.zmin(),
                        bbox.xmax(), bbox.ymax(), bbox.zmax(), t);
  }

  ///////////////////////////////////////////////////////////////////////////////
  /// Test if a particular point is outside of the object (Tree), knowing there is
  /// no intersection between its voxel and the tree.
  template<typename Mesh, typename Tree>
  bool is_outside_knowing_no_intersect(const internal::Hexmeshing::Point& p, const Tree& t)
  {
    CGAL::Side_of_triangle_mesh<Mesh, typename Kernel_traits
                                <internal::Hexmeshing::Point>::Kernel> s(t);
    CGAL::Bounded_side res=s(p);
    return res!=CGAL::ON_BOUNDED_SIDE; // && !=CGAL::ON_BOUNDARY ?
  }
}

#endif
