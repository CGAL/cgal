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

#ifndef CGAL_APOLLONIUS_GRAPH_2_ORIENTATION_NEW_2_H
#define CGAL_APOLLONIUS_GRAPH_2_ORIENTATION_NEW_2_H


// FIXME: We include the old traits class file for now to get the functors.
#include <CGAL/Apollonius_graph_traits_2.h>

CGAL_BEGIN_NAMESPACE


template <class K, class MTag>
class AG2_Orientation_test_new_2 : public AG2_Orientation_test_2<K, MTag>
{
public:
    typedef K                    Kernel;
    typedef typename K::FT       FT;
    typedef typename K::Site_2   Site_2;
    typedef typename K::Point_2  Point_2;

    typedef Orientation          result_type;
    typedef Arity_tag<3>         Arity;
    typedef Site_2               argument_type;

    Orientation operator() (const Site_2 &s0, const Site_2 &s1, const Site_2 &s2,
                            const Point_2 &q) const
    {
        FT x1 = s1.x() - s0.x(), y1 = s1.y() - s0.y(), w1 = s1.weight() - s0.weight();
        FT x2 = s2.x() - s0.x(), y2 = s2.y() - s0.y(), w2 = s2.weight() - s0.weight();
        FT xq =  q.x() - s0.x(), yq =  q.y() - s0.y();

        FT a1 = x1 * x1 + y1 * y1 - w1 * w1;
        FT a2 = x2 * x2 + y2 * y2 - w2 * w2;
        
        CGAL_assertion (sign (a1) > 0);
        CGAL_assertion (sign (a2) > 0);

        FT x = a1 * x2 - a2 * x1;
        FT y = a1 * y2 - a2 * y1;
        FT w = a1 * w2 - a2 * w1;
        FT s = x * xq + y * yq;

        Sign W = CGAL::sign (w);
        Sign S = CGAL::sign (s);

        if (W == 0) return Sign (- S);

        FT o = x * yq - y * xq;
        
        Sign O = CGAL::sign (o);

        if (S == 0) return Sign (O * W);

        if (W * S * O <= 0) return Sign (- S);

        FT i = w * w * (xq * xq + yq * yq) - s * s;

        Sign I = CGAL::sign (i);
        
        return Sign (S * I);
    }

};


CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_2_ORIENTATION_NEW_2_H
