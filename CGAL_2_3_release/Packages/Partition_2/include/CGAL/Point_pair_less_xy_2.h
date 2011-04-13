// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : include/CGAL/Point_pair_less_xy_2.h
// package       : $CGAL_Package: Partition_2 $
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// chapter       : Planar Polygon Partitioning
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: Lexicographical ordering of a pair of points
// ============================================================================

#ifndef CGAL_POINT_PAIR_LESS_XY_2_H
#define CGAL_POINT_PAIR_LESS_XY_2_H

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
     Point_pair_less_xy_2() : _less_xy_2(Traits().less_xy_2_object())
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

