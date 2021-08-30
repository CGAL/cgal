// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_TRAITS_WRAPPER_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_TRAITS_WRAPPER_2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>

namespace CGAL {

template<class Gt_base>
class Segment_Delaunay_graph_traits_wrapper_2 : public Gt_base
{
private:
  typedef Segment_Delaunay_graph_traits_wrapper_2<Gt_base> Self;

public:
  struct Triangle_2 {};

  Segment_Delaunay_graph_traits_wrapper_2() {}

  Segment_Delaunay_graph_traits_wrapper_2(const Gt_base&) {}
};

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_TRAITS_WRAPPER_2_H
