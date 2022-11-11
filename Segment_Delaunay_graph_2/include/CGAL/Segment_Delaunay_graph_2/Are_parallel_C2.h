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

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARE_PARALLEL_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARE_PARALLEL_C2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>
#include <CGAL/determinant.h>

namespace CGAL {

namespace SegmentDelaunayGraph_2 {

//-----------------------------------------------------------------------
//                           are parallel
//-----------------------------------------------------------------------

template< class K >
class Are_parallel_C2
{

public:
  typedef typename K::Site_2       Site_2;
  typedef typename K::Boolean      Boolean;
  typedef Boolean                  result_type;
  typedef Site_2                   argument_type;

private:
  typedef typename K::Segment_2    Segment_2;
  typedef typename K::FT           FT;

private:
  Boolean   predicate(const Site_2& p, const Site_2& q) const {
    CGAL_precondition( p.is_segment() && q.is_segment() );

    Segment_2 s1 = p.segment();
    Segment_2 s2 = q.segment();

    FT x1 = s1.source().x(),
      y1 = s1.source().y(),
      x2 = s1.target().x(),
      y2 = s1.target().y(),
      x3 = s2.source().x(),
      y3 = s2.source().y(),
      x4 = s2.target().x(),
      y4 = s2.target().y();

    FT det = determinant<FT>(x2 - x1, x4 - x3,
                             y2 - y1, y4 - y3);

    return ( CGAL::sign(det) == CGAL::ZERO );
  }

public:
  result_type operator()(const Site_2& p, const Site_2& q) const
  {
    return predicate(p, q);
  }
};

} //namespace SegmentDelaunayGraph_2

} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARE_PARALLEL_C2_H
