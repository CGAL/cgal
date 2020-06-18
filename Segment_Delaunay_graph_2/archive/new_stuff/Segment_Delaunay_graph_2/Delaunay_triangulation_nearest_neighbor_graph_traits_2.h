// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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

#ifndef CGAL_DELAUNAY_TRIANGULATION_NEAREST_NEIGHBOR_GRAPH_TRAITS_2_H
#define CGAL_DELAUNAY_TRIANGULATION_NEAREST_NEIGHBOR_GRAPH_TRAITS_2_H 1

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Voronoi_diagram_2/Site_accessors.h>

namespace CGAL {

template<class DG>
class Delaunay_triangulation_nearest_neighbor_graph_traits_2
{
public:
  typedef DG   Delaunay_graph;
  typedef typename Delaunay_graph::Geom_traits    Geom_traits;
  typedef typename Geom_traits::Point_2           Site_2;
  typedef typename Delaunay_graph::Vertex_handle  Vertex_handle;

  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Point_accessor<Site_2,Delaunay_graph,Tag_true>
  Access_site_2;

  Delaunay_triangulation_nearest_neighbor_graph_traits_2
  (const Geom_traits& gt = Geom_traits()) : gt_(gt) {}

#if 0
  struct Distance_comparator_2
  {
    //    typedef P                        Site_2;
    typedef CGAL::Comparison_result  result_type;

    result_type operator()(const Site_2& q, const Site_2& p1,
                           const Site_2& p2) const
    {
      // return SMALLER if p1 is closer to q than p2, EQUAL if p1 and p2
      // are at equal distance from q and LARGER is p2 is closer to q
      // than p1.
      typename Geom_traits::Compare_distance_2 comparator =
        Geom_traits().compare_distance_2_object();

      return comparator(q, p1, p2);
    }
  };
#else
  typedef typename Geom_traits::Compare_distance_2
  Distance_comparator_2;
#endif

  Access_site_2 access_site_2_object() const {
    return Access_site_2();
  }

  Distance_comparator_2 distance_comparator_2_object() const {
    return gt_.compare_distance_2_object();
    //Distance_comparator_2();
  }

private:
  const Geom_traits& gt_;
};

} //namespace CGAL

#endif // CGAL_DELAUNAY_TRIANGULATION_NEAREST_NEIGHBOR_GRAPH_TRAITS_2_H
