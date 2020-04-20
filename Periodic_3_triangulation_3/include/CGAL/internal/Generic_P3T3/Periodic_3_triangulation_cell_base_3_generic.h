// Copyright (c) 2020 XXXXX
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé
//                 Georg Osang

#ifndef CGAL_PERIODIC_3_TRIANGULATION_CELL_BASE_3_GENERIC_H
#define CGAL_PERIODIC_3_TRIANGULATION_CELL_BASE_3_GENERIC_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/Triangulation_cell_base_3.h>

#include <CGAL/Periodic_3_offset_3.h>

#include <array>
#include <fstream>

namespace CGAL {

template < typename Gt,
           typename Cb = Triangulation_cell_base_3<Gt> >
class Periodic_3_triangulation_cell_base_3_generic
  : public Cb
{
  typedef Cb                                            Base;
  typedef typename Base::Triangulation_data_structure   Tds;

public:
  typedef Gt                                            Geom_traits;
  typedef Tds                                           Triangulation_data_structure;
  typedef typename Tds::Vertex_handle                   Vertex_handle;
  typedef typename Tds::Cell_handle                     Cell_handle;
  typedef Periodic_3_offset_3                           Offset;

  template < typename TDS2 >
  struct Rebind_TDS
  {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other          Cb2;
    typedef Periodic_3_triangulation_cell_base_3_generic<Gt, Cb2>  Other;
  };

public:
  Periodic_3_triangulation_cell_base_3_generic() : Cb(), is_canonical(false) { }

  Periodic_3_triangulation_cell_base_3_generic(Vertex_handle v0,
                                               Vertex_handle v1,
                                               Vertex_handle v2,
                                               Vertex_handle v3)
    : Cb(v0, v1, v2, v3), is_canonical(false)
  { }

  Periodic_3_triangulation_cell_base_3_generic(Vertex_handle v0,
                                               Vertex_handle v1,
                                               Vertex_handle v2,
                                               Vertex_handle v3,
                                               Cell_handle n0,
                                               Cell_handle n1,
                                               Cell_handle n2,
                                               Cell_handle n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3), is_canonical(false)
  { }

  /// Periodic functions
  Offset offset(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i < 4 );
    return _off[i];
  }

  bool has_zero_offsets() const
  {
    return (_off[0] == Offset(0,0,0) &&
            _off[1] == Offset(0,0,0) &&
            _off[2] == Offset(0,0,0) &&
            _off[3] == Offset(0,0,0));
  }

  void set_offsets(const Offset& o0, const Offset& o1, const Offset& o2, const Offset& o3)
  {
    _off[0] = o0;
    _off[1] = o1;
    _off[2] = o2;
    _off[3] = o3;
  }

  void set_canonical_flag(const bool b) { is_canonical = b; }
  bool get_canonical_flag() const { return is_canonical; }

private:
  bool is_canonical;
  std::array<Offset, 4> _off;
};

template < class Tds >
inline
std::istream&
operator>>(std::istream &is, Periodic_3_triangulation_cell_base_3_generic<Tds>& )
// non combinatorial information. Default = nothing
{
  return is;
}

template < class Tds >
inline
std::ostream&
operator<<(std::ostream &os, const Periodic_3_triangulation_cell_base_3_generic<Tds>& )
// non combinatorial information. Default = nothing
{
  return os;
}

} //namespace CGAL

#endif //CGAL_PERIODIC_3_TRIANGULATION_CELL_BASE_3_GENERIC_H
