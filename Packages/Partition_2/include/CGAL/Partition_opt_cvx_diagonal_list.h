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
// file          : include/CGAL/Partition_opt_cvx_diagonal_list.h
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
// implementation: List of diagonals used in optimal convex partitioning alg.
// ============================================================================

#ifndef   CGAL_PARTITION_OPT_CVX_DIAGONAL_LIST_H
#define   CGAL_PARTITION_OPT_CVX_DIAGONAL_LIST_H

#include <utility>
#include <list>
#include <iostream>

typedef std::pair<int, int>                   Partition_opt_cvx_diagonal;
typedef std::list<Partition_opt_cvx_diagonal> Partition_opt_cvx_diagonal_list;

std::ostream& operator<<(std::ostream& os,
                         const Partition_opt_cvx_diagonal_list& d)
{
   Partition_opt_cvx_diagonal_list::const_iterator it;
   for (it = d.begin(); it != d.end(); it++)
   {
      os << "(" << (*it).first << ", " << (*it).second << ") ";
   }
   return os;
}

#endif // CGAL_PARTITION_OPT_CVX_DIAGONAL_LIST_H
