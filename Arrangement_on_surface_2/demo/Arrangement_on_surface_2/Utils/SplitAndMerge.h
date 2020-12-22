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

#ifndef ARRANGEMENT_DEMO_SPLIT_AND_MERGE
#define ARRANGEMENT_DEMO_SPLIT_AND_MERGE

// provides the same functionality as Arrangement_with_history_2::split_edge,
// but precompiled for all arrangements for better compilation speeds elsewhere
template <typename Arr_>
class Split_edge
{
public:
  using Arrangement = Arr_;
  using Traits = typename Arrangement::Geometry_traits_2;
  using Halfedge_handle = typename Arrangement::Halfedge_handle;
  using Point_2 = typename Traits::Point_2;
  using Equal_2 = typename Traits::Equal_2;

  void operator()(Arrangement*, Halfedge_handle, const Point_2&);
};

// provides the same functionality as Arrangement_with_history_2::merge_edge,
// but precompiled for all arrangements for better compilation speeds elsewhere
template <typename Arr_>
class Merge_edge
{
public:
  using Arrangement = Arr_;
  using Traits = typename Arrangement::Geometry_traits_2;
  using Halfedge_handle = typename Arrangement::Halfedge_handle;
  using Point_2 = typename Traits::Point_2;

  void mergeEdge(Arrangement*, Halfedge_handle, Halfedge_handle);
  bool areMergeable(Arrangement*, Halfedge_handle, Halfedge_handle);
};

#endif
