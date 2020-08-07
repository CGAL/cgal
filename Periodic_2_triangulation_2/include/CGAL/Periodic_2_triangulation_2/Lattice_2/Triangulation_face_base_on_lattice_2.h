// Copyright (c) 2020 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_P2T2_TRIANGULATION_FACE_BASE_ON_LATTICE_2_H
#define CGAL_P2T2_TRIANGULATION_FACE_BASE_ON_LATTICE_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_face_base_2.h>

#include <CGAL/Periodic_2_offset_2.h>

#include <array>
#include <fstream>

namespace CGAL {

template <typename Gt,
          typename Fb = Triangulation_face_base_2<Gt> >
class Triangulation_face_base_on_lattice_2
  : public Fb
{
  typedef Fb                                            Base;

public:
  typedef Gt                                            Geom_traits;
  typedef typename Fb::Vertex_handle                    Vertex_handle;
  typedef typename Fb::Face_handle                      Face_handle;
  typedef Periodic_2_offset_2                           Offset;

  template <typename TDS2>
  struct Rebind_TDS
  {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Triangulation_face_base_on_lattice_2<Gt, Fb2>  Other;
  };

public:
  Triangulation_face_base_on_lattice_2() : Fb(), is_canonical(false) { }

  Triangulation_face_base_on_lattice_2(Vertex_handle v0,
                                       Vertex_handle v1,
                                       Vertex_handle v2)
    : Fb(v0, v1, v2), is_canonical(false)
  { }

  Triangulation_face_base_on_lattice_2(Vertex_handle v0,
                                       Vertex_handle v1,
                                       Vertex_handle v2,
                                       Face_handle n0,
                                       Face_handle n1,
                                       Face_handle n2)
    : Fb(v0, v1, v2, n0, n1, n2), is_canonical(false)
  { }

  /// Periodic functions
  Offset offset(const int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i < 3 );
    return _off[i];
  }

  bool has_zero_offsets() const
  {
    return (_off[0] == Offset(0,0) &&
            _off[1] == Offset(0,0) &&
            _off[2] == Offset(0,0));
  }

  void set_offsets(const Offset& o0, const Offset& o1, const Offset& o2)
  {
    _off[0] = o0;
    _off[1] = o1;
    _off[2] = o2;
  }

  void set_canonical_flag(const bool b) { is_canonical = b; }
  bool get_canonical_flag() const { return is_canonical; }

private:
  bool is_canonical;
  std::array<Offset, 3> _off;
};

template <typename Gt, typename Fb>
inline std::istream& operator>>(std::istream &is, Triangulation_face_base_on_lattice_2<Gt, Fb>&)
{
  // non combinatorial information. Default = nothing
  return is;
}

template <typename Gt, typename Fb>
inline std::ostream& operator<<(std::ostream &os, const Triangulation_face_base_on_lattice_2<Gt, Fb>&)
{
  // non combinatorial information. Default = nothing
  return os;
}

} //namespace CGAL

#endif //CGAL_P2T2_TRIANGULATION_FACE_BASE_ON_LATTICE_2_H
