// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_hierarchy_vertex_base_2.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Triangulation_2
// author(s)     : Olivier Devillers <Olivivier.Devillers@sophia.inria.fr>
//                 Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_HIERARCHY_VERTEX_BASE_2_H
#define CGAL_TRIANGULATION_HIERARCHY_VERTEX_BASE_2_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

template < class Vbb>
class Triangulation_hierarchy_vertex_base_2
 : public Vbb
{
  typedef Vbb                                           Base;
  typedef typename Base::Triangulation_data_structure   Tds;

public:
  typedef typename Base::Point            Point;
  typedef Tds                             Triangulation_data_structure;
  typedef typename Tds::Vertex_handle     Vertex_handle;
  typedef typename Tds::Face_handle       Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vbb::template Rebind_TDS<TDS2>::Other      Vb2;
    typedef Triangulation_hierarchy_vertex_base_2<Vb2>          Other;
  };

  Triangulation_hierarchy_vertex_base_2()
    : Base(), _up(0), _down(0)
    {}
  Triangulation_hierarchy_vertex_base_2(const Point & p, Face_handle f)
    : Base(p,f), _up(0), _down(0)
    {}
  Triangulation_hierarchy_vertex_base_2(const Point & p)
    : Base(p), _up(0), _down(0)
    {}

  Vertex_handle up() {return _up;}
  Vertex_handle down() {return _down;}
  void set_up(Vertex_handle u) {_up=u;}
  void set_down(Vertex_handle d) {if (this) _down=d;}


 private:
  Vertex_handle  _up;    // same vertex one level above
  Vertex_handle  _down;  // same vertex one level below
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_HIERARCHY_VERTEX_BASE_2_H
