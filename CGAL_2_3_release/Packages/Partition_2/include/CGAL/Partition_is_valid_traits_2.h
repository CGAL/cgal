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
// file          : include/CGAL/Partition_is_valid_traits_2.h
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
// implementation: Partition validity traits class
// ============================================================================

#ifndef CGAL_PARTITION_2_IS_VALID_TRAITS_H
#define CGAL_PARTITION_2_IS_VALID_TRAITS_H

namespace CGAL {

template <class Traits, class PolygonIsValid>
class Partition_is_valid_traits_2 : public Traits
{
public:
   typedef typename Traits::Point_2         Point_2;
   typedef typename Traits::Polygon_2       Polygon_2;
   typedef typename Traits::Less_xy_2       Less_xy_2;
   typedef typename Traits::Leftturn_2      Leftturn_2;
   typedef typename Traits::Orientation_2   Orientation_2;

   typedef PolygonIsValid                   Is_valid;
   
   Is_valid
   is_valid_object(const Traits& traits) const
   {  return Is_valid(traits); }
};

}

#endif // CGAL_PARTITION_2_IS_VALID_TRAITS_H

