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
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_TRIANGULATION_FACE_BASE_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>


#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Dummy_tds_2.h>

namespace CGAL
{

template < typename Gt, typename Fb = Triangulation_face_base_2<Gt> >
class Periodic_2_triangulation_face_base_2
  : public Fb
{
  typedef Fb                                            Base;
  typedef typename Base::Triangulation_data_structure   Tds;

public:
  typedef Gt                                            Geom_traits;
  typedef Tds                                           Triangulation_data_structure;
  typedef typename Tds::Vertex_handle                   Vertex_handle;
  typedef typename Tds::Face_handle                     Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS
  {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Periodic_2_triangulation_face_base_2<Gt, Fb2>  Other;
  };

public:
  Periodic_2_triangulation_face_base_2()
    : Fb(), _off(0) {}

  Periodic_2_triangulation_face_base_2(Vertex_handle v0,
                                       Vertex_handle v1,
                                       Vertex_handle v2)
    : Fb(v0, v1, v2) , _off(0) {}

  Periodic_2_triangulation_face_base_2(Vertex_handle v0,
                                       Vertex_handle v1,
                                       Vertex_handle v2,
                                       Face_handle n0,
                                       Face_handle n1,
                                       Face_handle n2)
    : Fb(v0, v1, v2, n0, n1, n2), _off(0) {}

  /// Periodic functions
  int offset(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i < 3 );
    return ((_off >> 2 * i) & 3);
  }
  bool has_zero_offsets() const
  {
    return (_off & 63) == 0;
  }

  void set_offsets(unsigned int o0, unsigned int o1, unsigned int o2)
  {
    // 192=11000000
    _off = _off | 192;
    unsigned int off0[2] = {(o0 >> 1) & 1, (o0 & 1)};
    unsigned int off1[2] = {(o1 >> 1) & 1, (o1 & 1)};
    unsigned int off2[2] = {(o2 >> 1) & 1, (o2 & 1)};
    for (int i = 0; i < 2; i++)
      {
        unsigned int _off0 = ( _off    & 3);
        unsigned int _off1 = ((_off >> 2) & 3);
        unsigned int _off2 = ((_off >> 4) & 3);

        _off0 = ( (_off0 << 1) + off0[i]);
        _off1 = ( (_off1 << 1) + off1[i]);
        _off2 = ( (_off2 << 1) + off2[i]);

        // 252=11111100
        // 243=11110011
        // 207=11001111
        _off = ((_off & 252) | (_off0   ));
        _off = ((_off & 243) | (_off1 << 2));
        _off = ((_off & 207) | (_off2 << 4));
      }
  }

  void set_additional_flag(unsigned char b)
  {
    CGAL_assertion(b < 4);
    // 63=00111111
    _off = ((_off & 63) | (b << 6));
  }
  unsigned char get_additional_flag()
  {
    return (_off >> 6);
  }

private:
  // 2 respective bits are the _offset in x and y
  // right to left:
  // bit[0]-bit[1]: vertex(0),
  // bit[2]-bit[3]: vertex(1) and
  // bit[4]-bit[5]: vertex(2)
  // Thus the underlying data type needs to have at least 6 bit,
  // which is true for an unsigned char.
  // bit[6]: Used to convert 9 sheeted covering to a 1 sheeted covering
  unsigned char _off;
};

template < class Tds >
inline
std::istream&
operator>>(std::istream &is, Periodic_2_triangulation_face_base_2<Tds> &)
// non combinatorial information. Default = nothing
{
  return is;
}

template < class Tds >
inline
std::ostream&
operator<<(std::ostream &os, const Periodic_2_triangulation_face_base_2<Tds> &)
// non combinatorial information. Default = nothing
{
  return os;
}

} //namespace CGAL

#endif //CGAL_PERIODIC_2_TRIANGULATION_FACE_BASE_2_H
