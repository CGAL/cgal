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
// file          : include/CGAL/Indirect_not_less_yx_2.h
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
// implementation: Lexicographic comparison of points represented by iterators
// ============================================================================

#ifndef CGAL_INDIRECT_NOT_LESS_YX_2_H
#define CGAL_INDIRECT_NOT_LESS_YX_2_H

#include <CGAL/ch_utils.h>

namespace CGAL {

template <class Traits>
class Indirect_not_less_yx_2 
{
   public:
     typedef typename Traits::Less_yx_2     Less_yx_2;

     Indirect_not_less_yx_2(const Traits& traits) : 
         less_yx_2(traits.less_yx_2_object()) {}

     template <class Iterator>
     bool 
     operator()( const Iterator& p, const Iterator& q) const
     { return less_yx_2( *q, *p); }

   private:
     Less_yx_2 less_yx_2;
};

}

#endif // CGAL_INDIRECT_NOT_LESS_YX_2_H

