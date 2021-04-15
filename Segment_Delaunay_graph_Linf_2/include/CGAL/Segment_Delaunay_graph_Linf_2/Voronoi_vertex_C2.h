// Copyright (c) 2015  Università della Svizzera italiana.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VORONOI_VERTEX_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VORONOI_VERTEX_C2_H

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>



#include <CGAL/Segment_Delaunay_graph_Linf_2/basic.h>

#include <CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_ring_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_sqrt_field_new_C2.h>

namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

namespace Internal {

  template<class K,class M> struct Which_Voronoi_vertex_base_C2;

  template<class K>
  struct Which_Voronoi_vertex_base_C2<K,Integral_domain_without_division_tag>
  {
    typedef Voronoi_vertex_ring_C2<K>          Base;
  };

  template<class K>
  struct Which_Voronoi_vertex_base_C2<K,Field_with_sqrt_tag>
  {
    typedef Voronoi_vertex_sqrt_field_new_C2<K>    Base;
  };
} // namespace Internal

//----------------------------------------------------------------------

template<class K, class M>
class Voronoi_vertex_C2
  : public Internal::Which_Voronoi_vertex_base_C2<K,M>::Base
{
private:
  typedef typename Internal::Which_Voronoi_vertex_base_C2<K,M>::Base  Base;

protected:
  typedef typename Base::Site_2   Site_2;
public:
  Voronoi_vertex_C2(const Site_2& p, const Site_2& q,
                    const Site_2& r)
    : Base(p, q, r) {}
};


} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL



#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VORONOI_VERTEX_C2_H
