// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
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



#ifndef CGAL_APOLLONIUS_GRAPH_2_CONSTRUCTIONS_RTH2_H
#define CGAL_APOLLONIUS_GRAPH_2_CONSTRUCTIONS_RTH2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/determinant.h>

#include <CGAL/Apollonius_graph_2/Constructions_ftC2.h>

namespace CGAL {

template < class RT, class We >
inline
void
ad_circumcenterH2(const RT &hx1, const RT &hy1, const RT &hw1,
                  const We &w1,
                  const RT &hx2, const RT &hy2, const RT &hw2,
                  const We &w2,
                  const RT &hx3, const RT &hy3, const RT &hw3,
                  const We &w3,
                  RT       & cx, RT       & cy, RT       & cw)
{
  We cwt;
  ad_circumcircleH2(hx1, hy1, hw1, w1,
                    hx2, hy2, hw2, w2,
                    hx3, hy3, hw3, w3,
                    cx, cy, cw, cwt);
}

template < class RT, class We >
inline
void
ad_circumcircleH2(const RT &hx1, const RT &hy1, const RT &hw1,
                  const We &w1,
                  const RT &hx2, const RT &hy2, const RT &hw2,
                  const We &w2,
                  const RT &hx3, const RT &hy3, const RT &hw3,
                  const We &w3,
                  RT       & cx, RT       & cy, RT       & cw,
                  We      &cwt)
{
  ad_circumcircleC2<RT, We>(RT(hx1/hw1), RT(hy1/hw1), w1,
                            RT(hx2/hw2), RT(hy2/hw2), w2,
                            RT(hx3/hw3), RT(hy3/hw3), w3,
                            cx, cy, cwt);
  cw = RT(1);
}


template < class RT, class We >
void
ad_left_bitangent_lineH2(const FT &hx1, const FT &hy1, const We &w1,
                         const FT &hx2, const FT &hy2, const We &w2,
                         RT       &  a, RT       &  b, RT       & c)
{
  ad_left_bitangent_lineC2(RT(hx1/hw1), RT(hy1/hw1), w1,
                           RT(hx2/hw2), RT(hy2/hw2), w2,
                                     a,           b,  c);
}


} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_CONSTRUCTIONS_RTH2_H
