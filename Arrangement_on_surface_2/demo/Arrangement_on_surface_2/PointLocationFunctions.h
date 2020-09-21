// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ahmed Essam <theartful.ae@gmail.com>

#ifndef ARRANGEMENT_DEMO_POINT_LOCATION_FUNCTIONS
#define ARRANGEMENT_DEMO_POINT_LOCATION_FUNCTIONS

#include <CGAL/Object.h>
#include "ArrTraitsAdaptor.h"

class QPointF;

template <typename Arr_>
class PointLocationFunctions
{
public:
  using Arrangement = Arr_;
  using Traits = typename Arrangement::Geometry_traits_2;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Halfedge_around_vertex_const_circulator =
    typename Arrangement::Halfedge_around_vertex_const_circulator;
  using Kernel = typename ArrTraitsAdaptor<Traits>::Kernel;
  using Kernel_point_2 = typename Kernel::Point_2;

  // the QPointF versions are for convenience
  CGAL::Object locate(const Arrangement*, const Kernel_point_2&);
  CGAL::Object locate(const Arrangement*, const QPointF&);
  Face_const_handle getFace(const Arrangement*, const Kernel_point_2&);
  Face_const_handle getFace(const Arrangement*, const QPointF&);
  CGAL::Object rayShootUp(const Arrangement*, const QPointF&);
  CGAL::Object rayShootDown(const Arrangement*, const QPointF&);
};

#endif
