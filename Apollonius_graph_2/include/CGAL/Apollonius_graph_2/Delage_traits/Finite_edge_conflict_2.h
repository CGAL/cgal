// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>
//                 Christophe Delage <Christophe.Delage@sophia.inria.fr>
//                 David Millman <dlm336@cs.nyu.edu>

#ifndef CGAL_APOLLONIUS_GRAPH_2_FINITE_EDGE_CONFLICT_2_H
#define CGAL_APOLLONIUS_GRAPH_2_FINITE_EDGE_CONFLICT_2_H

#include <CGAL/Apollonius_graph_2/new_traits/Edge_conflict_2.h>

CGAL_BEGIN_NAMESPACE

CGAL_APOLLONIUS_GRAPH_2_BEGIN_NAMESPACE

//-----------------------------------------------------------------------
//                    Finite edge interior conflict
//-----------------------------------------------------------------------

template < class K, class Method_tag >
class Finite_edge_interior_conflict_new_2
  : public Edge_conflict_2<K, Method_tag>
{
private:
  typedef Edge_conflict_2<K,Method_tag>    Base;
public:
  typedef typename Base::Inverted_weighted_point   Inverted_weighted_point;
  typedef Weighted_point_inverter_2<K>             Weighted_point_inverter;
  typedef typename K::Site_2                       Site_2;
  typedef typename K::Point_2                      Point_2;
  typedef bool                                     result_type;

    inline
    bool operator()(const Site_2& p1, const Site_2& p2, 
		    const Site_2& q, bool b) const
    {
      Weighted_point_inverter inverter(p1);
      Point_2 origin(0,0);
      Site_2 origin_site(origin,0);
      return edge_conflict_test(inverter(p2),
				Inverted_weighted_point(origin_site,1),
				Inverted_weighted_point(origin_site,1),
				inverter(q), b, 2, 2);
    }

    inline
    bool operator()(const Site_2& p1, const Site_2& p2, 
		    const Site_2& p3, const Site_2& q, bool b) const
    {
      Weighted_point_inverter inverter(p2);
      Point_2 origin(0,0);
      Site_2 origin_site(origin,0);
      return edge_conflict_test(inverter(p1), 
				Inverted_weighted_point(origin_site,1),
				inverter(p3), inverter(q), b, 2, 1);
    }

    inline
    bool operator()(const Site_2& p1, const Site_2& p2,	const Site_2& p3,
		    const Site_2& p4, const Site_2& q, bool b) const
    {
      Weighted_point_inverter inverter(p2);
      return edge_conflict_test(inverter(p1), inverter(p4), inverter(p3), 
				inverter(q), b, 1, 1);
    }
};

CGAL_APOLLONIUS_GRAPH_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_2_FINITE_EDGE_CONFLICT_2_H
