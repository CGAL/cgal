// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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



#ifndef CGAL_APOLLONIUS_GRAPH_2_IS_DEGENERATE_EDGE_C2_H
#define CGAL_APOLLONIUS_GRAPH_2_IS_DEGENERATE_EDGE_C2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Apollonius_graph_2/basic.h>

#include <CGAL/Apollonius_graph_2/Predicate_constructions_C2.h>

#include <CGAL/Apollonius_graph_2/Incircle_C2.h>
#include <CGAL/Apollonius_graph_2/Finite_edge_test_C2.h>

namespace CGAL {

namespace ApolloniusGraph_2 {


//--------------------------------------------------------------------

template < class K, class MTag >
class Is_degenerate_edge_2
{
public:
  typedef K                                 Kernel;
  typedef MTag                              Method_tag;

  typedef typename K::Site_2                Site_2;
  typedef Weighted_point_inverter_2<K>      Weighted_point_inverter;
  typedef Inverted_weighted_point_2<K>      Inverted_weighted_point;
  typedef Bitangent_line_2<K>               Bitangent_line;
  typedef Voronoi_circle_2<K>               Voronoi_circle;
  typedef typename K::FT                    FT;
  typedef typename K::Sign                  Sign;
  typedef typename K::Comparison_result     Comparison_result;

  typedef Order_on_finite_bisector_2<K>     Order_on_finite_bisector;

  typedef Sign_of_distance_from_CCW_circle_2<K>
                                          Sign_of_distance_from_CCW_circle;

public:
  typedef Site_2             argument_type;
  typedef bool               result_type;

  bool operator()(const Site_2& p1, const Site_2& p2,
                  const Site_2& p3, const Site_2& p4) const
  {
    Method_tag tag;

    Weighted_point_inverter inverter(p1);
    Inverted_weighted_point u2 = inverter(p2);
    Inverted_weighted_point u3 = inverter(p3);
    Inverted_weighted_point u4 = inverter(p4);

    Sign s;

    Bitangent_line blinv_23(u2, u3);
    s = Sign_of_distance_from_CCW_circle()(blinv_23, u4, tag);
    if ( s != ZERO ) { return false; }

    Bitangent_line blinv_42(u4, u2);
    s = Sign_of_distance_from_CCW_circle()(blinv_42, u3, tag);
    if ( s != ZERO ) { return false; }

    Voronoi_circle vc_123(blinv_23);
    Voronoi_circle vc_142(blinv_42);
    Comparison_result r =
      Order_on_finite_bisector()(vc_123, vc_142, p1, p2, tag);

    return ( r == EQUAL );
  }
};

//--------------------------------------------------------------------

} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_IS_DEGENERATE_EDGE_C2_H
