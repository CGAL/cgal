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
// file          : include/CGAL/Turn_reverser.h
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
// implementation: Function object to reverse a turn predicate
// ============================================================================

#ifndef   CGAL_TURN_REVERSER_H
#define   CGAL_TURN_REVERSER_H

namespace CGAL {

template <class Point_2, class TurnPredicate>
class Turn_reverser 
{
public:
   Turn_reverser() {}
   Turn_reverser( const TurnPredicate& t ): turn(t) {}

   bool operator() (const Point_2& p1, 
                    const Point_2& p2, 
                    const Point_2& p3) const
   {   return turn(p2, p1, p3); }

private:
   TurnPredicate turn;
};


}

#endif // CGAL_TURN_REVERSER_H
