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
#ifndef CGAL_HALFEDGEDS_ITEMS_H
#include <CGAL/HalfedgeDS_items.h>
#endif
#ifndef CGAL_HALFEDGEDS_USING_IN_PLACE_LIST_H
#include <CGAL/HalfedgeDS_using_in_place_list.h>
#endif

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
    template <class p_Traits, class p_Items = HalfedgeDS_items>
    class HalfedgeDS_default
        : public HalfedgeDS_using_in_place_list< p_Traits, p_Items> {
    public:
        typedef p_Traits Traits;
        typedef HalfedgeDS_using_in_place_list<p_Traits, p_Items> DS;
        typedef typename DS::size_type size_type;
        HalfedgeDS_default() {}
        HalfedgeDS_default( size_type v, size_type h, size_type f)
            : HalfedgeDS_using_in_place_list< p_Traits, p_Items>(v,h,f) {}
    };
    #define CGAL_HALFEDGEDS_DEFAULT  ::CGAL::HalfedgeDS_default
#else
    struct HalfedgeDS_default {
      template <class p_Traits, class p_Items = HalfedgeDS_items>
      class HDS
          : public HalfedgeDS_using_in_place_list::HDS<p_Traits, p_Items> {
      public:
          typedef p_Traits Traits;
          typedef HalfedgeDS_using_in_place_list::HDS<p_Traits, p_Items> DS;
          typedef typename DS::size_type size_type;
          //HDS() {}
          //HDS( size_type v, size_type h, size_type f)
          // : HalfedgeDS_using_in_place_list::HDS<p_Traits, p_Items>(v,h,f) {}
          HDS();
          HDS( size_type v, size_type h, size_type f);
      };
    };
    template <class p_Traits, class p_Items>
    HalfedgeDS_default::HDS<p_Traits, p_Items>:: HDS() {}

    template <class p_Traits, class p_Items>
    HalfedgeDS_default::HDS<p_Traits, p_Items>::
    HDS( size_type v, size_type h, size_type f)
        : HalfedgeDS_using_in_place_list::HDS<p_Traits, p_Items>(v,h,f) {}

    #define CGAL_HALFEDGEDS_DEFAULT  ::CGAL::HalfedgeDS_default::HDS
#endif // CGAL_CFG_NO_TMPL_IN_TMPL_PARAM //

    CGAL_END_NAMESPACE
#endif // CGAL_HALFEDGEDS_DEFAULT_H //
// EOF //
