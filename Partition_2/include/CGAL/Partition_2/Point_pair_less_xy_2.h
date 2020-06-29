// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_POINT_PAIR_LESS_XY_2_H
#define CGAL_POINT_PAIR_LESS_XY_2_H

#include <CGAL/license/Partition_2.h>


#include <utility>

namespace CGAL {

//
// Given two pairs of points determine their lexicographical order by first
// comparing the first points lexicographically and then the second points if
// the first are equal
//
template <class Traits>
class Point_pair_less_xy_2
{
   typedef typename Traits::Point_2           Point_2;
   typedef std::pair<Point_2, Point_2>        Point_pair;
   typedef typename Traits::Less_xy_2         Less_xy_2;

   public:
     Point_pair_less_xy_2(const Traits& traits) : _less_xy_2(traits.less_xy_2_object())
     { }


     bool
     operator()(const Point_pair& p, const Point_pair& q) const
     {
        if (_less_xy_2(p.first, q.first))
            return true;
        else if (_less_xy_2(q.first, p.first))
            return false;
        else if (_less_xy_2(p.second, q.second))
            return true;
        else
            return false;
     }

   private:
      Less_xy_2 _less_xy_2;
};

}

#endif // CGAL_POINT_PAIR_LESS_XY_H
