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
// package       : $CGAL_Package: HalfedgeDS 3.3 (27 Sep 2000) $
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

#include <CGAL/HalfedgeDS_items_2.h>
#include <CGAL/HalfedgeDS_list.h>
#include <CGAL/memory.h>

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM

template <class Traits_, class HalfedgeDSItems = HalfedgeDS_items_2, 
          class Alloc = CGAL_ALLOCATOR(int)>
class HalfedgeDS_default 
    : public HalfedgeDS_list< Traits_, HalfedgeDSItems, Alloc> {
public:
    typedef Traits_                                          Traits;
    typedef HalfedgeDS_list<Traits_, HalfedgeDSItems, Alloc> DS;
    typedef typename DS::size_type                           size_type;
    HalfedgeDS_default() {}
    HalfedgeDS_default( size_type v, size_type h, size_type f)
        : HalfedgeDS_list< Traits_, HalfedgeDSItems, Alloc>(v,h,f) {}
};
#define CGAL_HALFEDGEDS_DEFAULT  ::CGAL::HalfedgeDS_default

#else //  CGAL_CFG_NO_TMPL_IN_TMPL_PARAM //

struct HalfedgeDS_default {
  template <class Traits_, class HalfedgeDSItems = HalfedgeDS_items_2, 
            class Alloc = CGAL_ALLOCATOR(int)>
  class HDS : public HalfedgeDS_list::HDS<Traits_, HalfedgeDSItems, Alloc> {
  public:
      typedef Traits_                                               Traits;
      typedef HalfedgeDS_list::HDS<Traits_, HalfedgeDSItems, Alloc> DS;
      typedef typename DS::size_type                                size_type;
      HDS();
      HDS( size_type v, size_type h, size_type f);
  };
};

template <class Traits_, class HalfedgeDSItems, class Alloc>
HalfedgeDS_default::HDS<Traits_, HalfedgeDSItems, Alloc>:: HDS() {}

template <class Traits_, class HalfedgeDSItems, class Alloc>
HalfedgeDS_default::HDS<Traits_, HalfedgeDSItems, Alloc>::
HDS( size_type v, size_type h, size_type f)
    : HalfedgeDS_list::HDS<Traits_, HalfedgeDSItems, Alloc>(v,h,f) {}

#define CGAL_HALFEDGEDS_DEFAULT  ::CGAL::HalfedgeDS_default::HDS

#endif // CGAL_CFG_NO_TMPL_IN_TMPL_PARAM //

CGAL_END_NAMESPACE
#endif // CGAL_HALFEDGEDS_DEFAULT_H //
// EOF //
