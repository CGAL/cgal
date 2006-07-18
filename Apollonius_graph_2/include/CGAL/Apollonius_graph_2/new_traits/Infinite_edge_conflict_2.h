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

#ifndef CGAL_APOLLONIUS_GRAPH_2_INFINITE_EDGE_CONFLICT_2_H
#define CGAL_APOLLONIUS_GRAPH_2_INFINITE_EDGE_CONFLICT_2_H


// FIXME: We include the old traits class file for now to get the functors.
#include <CGAL/Apollonius_graph_traits_2.h>

#include <CGAL/Apollonius_graph_2/new_traits/Edge_conflict_2.h>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------
//                   Infinite edge interior conflict
//-----------------------------------------------------------------------

template < class K, class Method_tag >
class Infinite_edge_interior_conflict_new_2 : public Edge_conflict_2<K, Method_tag>
{
public:
    typedef CGAL::Weighted_point_inverter<K> Weighted_point_inverter;
    typedef CGAL::Inverted_weighted_point<K> Inverted_weighted_point;
    typedef typename K::Site_2               Site_2;
    typedef typename K::Point_2              Point_2;
    typedef bool                             result_type;
    typedef Arity_tag<5>                     Arity;

    inline
    bool operator()(const Site_2& p2, const Site_2& p3, 
            const Site_2& p4, const Site_2& q, bool b) const
    {
        Weighted_point_inverter inverter(p2);
        return edge_conflict_test(
                Inverted_weighted_point(Site_2(Point_2(0,0),0), 1),
                inverter(p4), inverter(p3), inverter(q), b, 1, 1);
    }
};


CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_2_INFINITE_EDGE_CONFLICT_2_H
