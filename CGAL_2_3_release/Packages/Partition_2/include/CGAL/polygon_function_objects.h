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
// file          : include/CGAL/polygon_function_objects.h
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
// implementation: Function objects for testing polygon properties
// ============================================================================

#ifndef POLYGON_FUNCTION_OBJECTS_H
#define POLYGON_FUNCTION_OBJECTS_H

#include<CGAL/is_y_monotone_2.h>

namespace CGAL {

template <class Traits>
class Is_vacuously_valid 
{
  public:

     Is_vacuously_valid(Traits ) {}

     template <class ForwardIterator>
     bool operator()(ForwardIterator, ForwardIterator)
     {  return true; }

};


template <class Traits>
class Is_convex_2
{
  public:
     Is_convex_2(Traits t): traits(t) {}
  
     template <class ForwardIterator>
     bool operator()(ForwardIterator first, ForwardIterator last)
     {  return is_convex_2(first, last, traits); }

  private:
     Traits  traits;
};

template <class Traits>
class Is_y_monotone_2
{
  public:
     Is_y_monotone_2(Traits t): traits(t) {}
  
     template <class ForwardIterator>
     bool operator()(ForwardIterator first, ForwardIterator last)
     {  return is_y_monotone_2(first, last, traits); }

  private:
     Traits  traits;
};

}

#endif // POLYGON_FUNCTION_OBJECTS_H
