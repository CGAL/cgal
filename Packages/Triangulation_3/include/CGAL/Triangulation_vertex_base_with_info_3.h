// ============================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_vertex_base_with_info_3.h
// revision      : $Revision$
// author(s)     : Sylvain Pion
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_VERTEX_BASE_WITH_INFO_3_H
#define CGAL_TRIANGULATION_VERTEX_BASE_WITH_INFO_3_H

#include <CGAL/Triangulation_vertex_base_3.h>

CGAL_BEGIN_NAMESPACE

template < typename Info_, typename GT,
           typename Vb = Triangulation_vertex_base_3<GT> >
class Triangulation_vertex_base_with_info_3
  : public Vb
{
  Info_ _info;
public:

  typedef typename Vb::Cell_handle                   Cell_handle;
  typedef typename Vb::Point                         Point;
  typedef Info_                                      Info;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other          Vb2;
    typedef Triangulation_vertex_base_with_info_3<Info, GT, Vb2>   Other;
  };

  Triangulation_vertex_base_with_info_3()
    : Vb() {}

  Triangulation_vertex_base_with_info_3(const Point & p)
    : Vb(p) {}

  Triangulation_vertex_base_with_info_3(const Point & p, Cell_handle c)
    : Vb(p, c) {}

  Triangulation_vertex_base_with_info_3(Cell_handle c)
    : Vb(c) {}

  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_VERTEX_BASE_WITH_INFO_3_H
