// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>
//                 Christophe Delage <Christophe.Delage@sophia.inria.fr>
//                 David Millman <dlm336@cs.nyu.edu>

#ifndef CGAL_APOLLONIUS_GRAPH_2_INFINITE_EDGE_CONFLICT_2_H
#define CGAL_APOLLONIUS_GRAPH_2_INFINITE_EDGE_CONFLICT_2_H


#include <CGAL/Apollonius_graph_2/Delage_traits/Edge_conflict_2.h>

namespace CGAL {

namespace ApolloniusGraph_2 {

//-----------------------------------------------------------------------
//                   Infinite edge interior conflict
//-----------------------------------------------------------------------

template < class K, class Method_tag >
class Infinite_edge_interior_conflict_new_2
  : public Edge_conflict_2<K, Method_tag>
{
private:
  typedef Edge_conflict_2<K,Method_tag>              Base;
  using Base::edge_conflict_test;

public:
  typedef Weighted_point_inverter_2<K>               Weighted_point_inverter;
  typedef typename Base::Inverted_weighted_point     Inverted_weighted_point;
  typedef typename K::Site_2                         Site_2;
  typedef typename K::Point_2                        Point_2;
  typedef bool                                       result_type;

  inline
  bool operator()(const Site_2& p2, const Site_2& p3, 
		    const Site_2& p4, const Site_2& q, bool b) const
  {
    Weighted_point_inverter inverter(p2);
    Point_2 origin(0,0);
    Site_2 origin_site(origin,0);
    return edge_conflict_test(Inverted_weighted_point(origin_site, 1),
			      inverter(p4), inverter(p3), inverter(q),
			      b, 1, 1);
  }
};


} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_INFINITE_EDGE_CONFLICT_2_H
