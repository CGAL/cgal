// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_vertex_base_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_VERTEX_BASE_3_H
#define CGAL_TRIANGULATION_VERTEX_BASE_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Triangulation_ds_vertex_base_3.h>

CGAL_BEGIN_NAMESPACE

template < typename GT, typename DSVb = Triangulation_ds_vertex_base_3<> >
class Triangulation_vertex_base_3
  : public DSVb
{
public:

  typedef typename DSVb::Cell_handle                   Cell_handle;
  typedef GT                                           Geom_traits;
  typedef typename GT::Point_3                         Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename DSVb::template Rebind_TDS<TDS2>::Other  DSVb2;
    typedef Triangulation_vertex_base_3<GT, DSVb2>           Other;
  };

  Triangulation_vertex_base_3()
    : DSVb() {}

  Triangulation_vertex_base_3(const Point & p)
    : DSVb(), _p(p) {}

  Triangulation_vertex_base_3(const Point & p, Cell_handle c)
    : DSVb(c), _p(p) {}

  Triangulation_vertex_base_3(Cell_handle c)
    : DSVb(c), _p() {}

  const Point & point() const
  { return _p; }

  Point & point()
  { return _p; }

  void set_point(const Point & p)
  { _p = p; }

private:
  Point _p;
};

template < class GT, class DSVb >
std::istream&
operator>>(std::istream &is, Triangulation_vertex_base_3<GT, DSVb> &v)
  // non combinatorial information. Default = point
{
  return is >> static_cast<DSVb&>(v) >> v.point();
}

template < class GT, class DSVb >
std::ostream&
operator<<(std::ostream &os, const Triangulation_vertex_base_3<GT, DSVb> &v)
  // non combinatorial information. Default = point
{
  return os << static_cast<const DSVb&>(v) << v.point();
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_VERTEX_BASE_3_H
