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

#include "SplitAndMerge.h"
#include "ArrangementTypes.h"

template <typename Arr_>
void Split_edge<Arr_>::operator()(
  Arrangement* arr, Halfedge_handle hei, const Point_2& pt)
{
  // check if split point lies on the edge of the curve
  Equal_2 areEqual(arr->traits()->equal_2_object());
  if (
    (!hei->source()->is_at_open_boundary() &&
     areEqual(hei->source()->point(), pt)) ||
    (!hei->target()->is_at_open_boundary() &&
     areEqual(hei->target()->point(), pt)))
  { return; }

  arr->split_edge(hei, pt);
}

template <typename Arr_>
void Merge_edge<Arr_>::mergeEdge(
  Arrangement* arr, Halfedge_handle h1, Halfedge_handle h2)
{
  arr->merge_edge(h1, h2);
}

template <typename Arr_>
bool Merge_edge<Arr_>::areMergeable(
  Arrangement* arr, Halfedge_handle h1, Halfedge_handle h2)
{
  return arr->are_mergeable(h1, h2);
}

ARRANGEMENT_DEMO_SPECIALIZE_ARR(Split_edge)
ARRANGEMENT_DEMO_SPECIALIZE_ARR(Merge_edge)
