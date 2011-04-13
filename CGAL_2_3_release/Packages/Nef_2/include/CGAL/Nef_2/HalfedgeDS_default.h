// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : HalfedgeDS_default.h
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: HalfedgeDS_2 3.1 (26 Mar 1999) $
// source        : hds.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Halfedge Data Structure Default Implementation for CGAL.
// ============================================================================

#ifndef CGAL_HALFEDGEDS_DEFAULT_H
#define CGAL_HALFEDGEDS_DEFAULT_H 1

#include <CGAL/Nef_2/HalfedgeDS_using_in_place_list.h>

CGAL_BEGIN_NAMESPACE

template <class p_Traits, class p_Items> // = HalfedgeDS_items>
class HalfedgeDS_default
  : public HalfedgeDS_using_in_place_list< p_Traits, p_Items> {
public:
  typedef p_Traits Traits;
  typedef size_t size_type;
  HalfedgeDS_default() {}
  HalfedgeDS_default( size_type v, size_type h, size_type f)
    : HalfedgeDS_using_in_place_list< p_Traits, p_Items>(v,h,f) {}
};

CGAL_END_NAMESPACE

#endif // CGAL_HALFEDGEDS_DEFAULT_H //
// EOF //
