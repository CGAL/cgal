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
// file          : include/CGAL/Triangulation_cell_base_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

// cell of a triangulation of any dimension <=3

#ifndef CGAL_TRIANGULATION_CELL_BASE_3_H
#define CGAL_TRIANGULATION_CELL_BASE_3_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Triangulation_ds_cell_base_3.h>

CGAL_BEGIN_NAMESPACE

template < typename GT, typename Cb = Triangulation_ds_cell_base_3<> >
class Triangulation_cell_base_3
  : public Cb
{
public:
  typedef typename Cb::Vertex_handle                   Vertex_handle;
  typedef typename Cb::Cell_handle                     Cell_handle;

  typedef GT                                           Geom_traits;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other  Cb2;
    typedef Triangulation_cell_base_3<GT, Cb2>             Other;
  };

  Triangulation_cell_base_3()
    : Cb() {}

  Triangulation_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                            Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {}

  Triangulation_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                            Vertex_handle v2, Vertex_handle v3,
                            Cell_handle   n0, Cell_handle   n1,
                            Cell_handle   n2, Cell_handle   n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3) {}
};

// The following should be useless.
#if 0
template < class GT, class Cb >
inline
std::istream&
operator>>(std::istream &is, Triangulation_cell_base_3<GT, Cb> &c)
  // non combinatorial information. Default = nothing
{
  return is >> static_cast<Cb&>(c);
}

template < class GT, class Cb >
inline
std::ostream&
operator<<(std::ostream &os, const Triangulation_cell_base_3<GT, Cb> &c)
  // non combinatorial information. Default = nothing
{
  return os << static_cast<const Cb&>(c);
}
#endif

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_CELL_BASE_3_H
