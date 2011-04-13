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
// file          : include/CGAL/Circulator_list.h
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
// implementation: List of circulators
// ============================================================================
#ifndef CGAL_CIRCULATOR_LIST_H
#define CGAL_CIRCULATOR_LIST_H

#include <list>
#include <iostream>

namespace CGAL {

template <class Circulator>
class Circulator_list : public std::list<Circulator> 
{

public:
  Circulator_list() 
  {}

  Circulator_list(Circulator first) 
  {
      if (first == NULL) return;

      Circulator current = first;
      do 
      {
         push_back(current);
      } while (++current != first);
  }
};

template <class C>
std::ostream& operator<<(std::ostream& os, const Circulator_list<C>& c)
{
     typename Circulator_list<C>::const_iterator current;
     for (current = c.begin(); current != c.end(); current++)
     {
        os << **current << " ";
     }
     return os;
}
}

#endif // CGAL_CIRCULATOR_LIST_H
