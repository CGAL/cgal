// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_VERTEX_BASE_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>


#include <CGAL/basic.h>
#include <CGAL/Dummy_tds_2.h>
#include <CGAL/Periodic_2_offset_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>

namespace CGAL
{


template < class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
class Periodic_2_triangulation_vertex_base_2
  : public Vb
{
  typedef Vb                                            Base;
  typedef typename Vb::Triangulation_data_structure     Tds;

public:
  typedef Gt                                            Geom_traits;
  typedef Tds                                           Triangulation_data_structure;

  typedef typename Tds::Vertex_handle                   Vertex_handle;
  typedef typename Tds::Face_handle                     Face_handle;
  typedef typename Gt::Point_2                          Point;
  typedef Periodic_2_offset_2                           Offset;

  template < typename Tds2 >
  struct Rebind_TDS
  {
    typedef typename Vb::template Rebind_TDS<Tds2>::Other      Vb2;
    typedef Periodic_2_triangulation_vertex_base_2<Gt, Vb2>    Other;
  };

public:
  Periodic_2_triangulation_vertex_base_2() : Base() {}
  Periodic_2_triangulation_vertex_base_2(const Point & p)
    : Base(p), _off(), _offset_flag(false) {}
  Periodic_2_triangulation_vertex_base_2(const Point & p, Face_handle f)
    : Base(f, p), _off(), _offset_flag(false) {}
  Periodic_2_triangulation_vertex_base_2(Face_handle f)
    : Base(f), _off(), _offset_flag(false) {}

  const Offset& offset() const
  {
    return _off;
  }

  void set_offset(const Offset& off)
  {
    _off = off;
    _offset_flag = true;
  }

  void clear_offset()
  {
    _offset_flag = false;
    _off = Offset();
  }

  bool get_offset_flag() const
  {
    return _offset_flag;
  }

private:
  /// The offset is needed to be able to copy a triangulation that is
  /// not on the 1-cover.

  /// Normal copying of the vertices would give multiple vertices with
  /// the same location (the periodic copies) and we wouldn't be able
  /// to distinguish them anymore.
  Offset _off;
  /// The flag is used to test whether _off has been set.
  bool _offset_flag;
};

template < class Tds >
inline
std::istream&
operator>>(std::istream &is, Periodic_2_triangulation_vertex_base_2<Tds> &)
// no combinatorial information.
{
  return is;
}

template < class Tds >
inline
std::ostream&
operator<<(std::ostream &os,
           const Periodic_2_triangulation_vertex_base_2<Tds> &)
// no combinatorial information.
{
  return os;
}

} //namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_VERTEX_BASE_2_H
