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
// file          : include/CGAL/Indirect_less_xy_2.h
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
// implementation: Given iterators to two points, compares them using Less_xy
// ============================================================================

#ifndef CGAL_INDIRECT_LESS_XY_2_H
#define CGAL_INDIRECT_LESS_XY_2_H

namespace CGAL {

template <class Traits>
class Indirect_less_xy_2 
{
   public:
     typedef typename Traits::Less_xy_2     Less_xy_2;

     Indirect_less_xy_2() : _less_xy_2(Traits().less_xy_2_object()) 
     { }
     
     template <class Iterator>
     bool 
     operator()(Iterator p, Iterator q) const
     { 
        return _less_xy_2(*p, *q);
     }

   private:
     Less_xy_2 _less_xy_2;
};

}

#endif // CGAL_INDIRECT_LESS_XY_2_H

