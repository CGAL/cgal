// Copyright (c) 1999-2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

// cell of a triangulation of any dimension <=3

#ifndef CGAL_PERIODIC_3_TRIANGULATION_DS_CELL_BASE_3_H
#define CGAL_PERIODIC_3_TRIANGULATION_DS_CELL_BASE_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/internal/Dummy_tds_3.h>

namespace CGAL {

template < typename TDS = void >
class Periodic_3_triangulation_ds_cell_base_3
{
public:
  typedef TDS                          Triangulation_data_structure;
  typedef typename TDS::Vertex_handle  Vertex_handle;
  typedef typename TDS::Cell_handle    Cell_handle;
  typedef typename TDS::Vertex         Vertex;
  typedef typename TDS::Cell           Cell;
  typedef typename TDS::Cell_data      TDS_data;

  template <typename TDS2>
  struct Rebind_TDS {
    typedef Periodic_3_triangulation_ds_cell_base_3<TDS2> Other;
  };

  Periodic_3_triangulation_ds_cell_base_3() : _additional_flag(0), off(0) {}

  Periodic_3_triangulation_ds_cell_base_3(
      const Vertex_handle& v0, const Vertex_handle& v1,
      const Vertex_handle& v2, const Vertex_handle& v3)
    : V{v0, v1, v2, v3},
      _additional_flag(0), off(0) {}

  Periodic_3_triangulation_ds_cell_base_3(
      const Vertex_handle& v0, const Vertex_handle& v1,
      const Vertex_handle& v2, const Vertex_handle& v3,
      const Cell_handle&   n0, const Cell_handle&   n1,
      const Cell_handle&   n2, const Cell_handle&   n3)
    : N{n0, n1, n2, n3},
      V{v0, v1, v2, v3},
      _additional_flag(0), off(0) {}

  // ACCESS FUNCTIONS

  const Vertex_handle& vertex(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 3 );
    CGAL_assume( i >= 0 && i <= 3 );
    return V[i];
  }

  bool has_vertex(const Vertex_handle& v) const
  {
    return (V[0] == v) || (V[1] == v) || (V[2]== v) || (V[3]== v);
  }

  bool has_vertex(const Vertex_handle& v, int & i) const
    {
      if (v == V[0]) { i = 0; return true; }
      if (v == V[1]) { i = 1; return true; }
      if (v == V[2]) { i = 2; return true; }
      if (v == V[3]) { i = 3; return true; }
      return false;
    }

  int index(const Vertex_handle& v) const
  {
    if (v == V[0]) { return 0; }
    if (v == V[1]) { return 1; }
    if (v == V[2]) { return 2; }
    CGAL_triangulation_assertion( v == V[3] );
    return 3;
  }

  const Cell_handle& neighbor(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 3);
    return N[i];
  }

  bool has_neighbor(const Cell_handle& n) const
  {
    return (N[0] == n) || (N[1] == n) || (N[2] == n) || (N[3] == n);
  }

  bool has_neighbor(const Cell_handle& n, int & i) const
  {
    if(n == N[0]){ i = 0; return true; }
    if(n == N[1]){ i = 1; return true; }
    if(n == N[2]){ i = 2; return true; }
    if(n == N[3]){ i = 3; return true; }
    return false;
  }

  int index(const Cell_handle& n) const
  {
    if (n == N[0]) return 0;
    if (n == N[1]) return 1;
    if (n == N[2]) return 2;
    CGAL_triangulation_assertion( n == N[3] );
    return 3;
  }

  int offset(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 3 );
    return ((off>>3*i)&7);
  }

  // SETTING

  void set_vertex(int i, const Vertex_handle& v)
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 3);
    V[i] = v;
  }

  void set_neighbor(int i, const Cell_handle& n)
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 3);
    N[i] = n;
  }

  void set_vertices()
  {
    V[0] = V[1] = V[2] = V[3] = Vertex_handle();
  }

  void set_vertices(const Vertex_handle& v0, const Vertex_handle& v1,
                    const Vertex_handle& v2, const Vertex_handle& v3)
  {
    V[0] = v0;
    V[1] = v1;
    V[2] = v2;
    V[3] = v3;
  }

  void set_neighbors()
  {
    N[0] = N[1] = N[2] = N[3] = Cell_handle();
  }

  void set_neighbors(const Cell_handle& n0, const Cell_handle& n1,
                     const Cell_handle& n2, const Cell_handle& n3)
  {
    N[0] = n0;
    N[1] = n1;
    N[2] = n2;
    N[3] = n3;
  }

  void set_offset(const Vertex_handle vh, int o)
  {
    CGAL_triangulation_precondition(has_vertex(vh));
    int vhi = index(vh);

    unsigned int offo[3] = {static_cast<unsigned int>((o&1)),
                            static_cast<unsigned int>((o>>1)&1),
                            static_cast<unsigned int>((o>>2)&1)};

    int bit_offset = 3 * vhi;

    // first reset the bit to 0 (AND), then assign the value given in input (OR)
    off = (off & ~(1 <<  bit_offset))      | (offo[0] <<  bit_offset);
    off = (off & ~(1 << (bit_offset + 1))) | (offo[1] << (bit_offset + 1));
    off = (off & ~(1 << (bit_offset + 2))) | (offo[2] << (bit_offset + 2));

    CGAL_postcondition(offset(vhi) == o);
  }

  void set_offsets(int o0,int o1,int o2,int o3)
  {
    off = 0;
    // The following explicit cast are needed according to the Intel
    // Compiler version 12.
    unsigned int off0[3] = {static_cast<unsigned int>((o0>>2)&1),
                            static_cast<unsigned int>((o0>>1)&1),
                            static_cast<unsigned int>((o0&1))};
    unsigned int off1[3] = {static_cast<unsigned int>((o1>>2)&1),
                            static_cast<unsigned int>((o1>>1)&1),
                            static_cast<unsigned int>((o1&1))};
    unsigned int off2[3] = {static_cast<unsigned int>((o2>>2)&1),
                            static_cast<unsigned int>((o2>>1)&1),
                            static_cast<unsigned int>((o2&1))};
    unsigned int off3[3] = {static_cast<unsigned int>((o3>>2)&1),
                            static_cast<unsigned int>((o3>>1)&1),
                            static_cast<unsigned int>((o3&1))};
    for (int i=0; i<3; i++) {
      unsigned int _off0 = (off&7);
      unsigned int _off1 = ((off>>3)&7);
      unsigned int _off2 = ((off>>6)&7);
      unsigned int _off3 = ((off>>9)&7);

      _off0 = ( (_off0<<1) + off0[i]);
      _off1 = ( (_off1<<1) + off1[i]);
      _off2 = ( (_off2<<1) + off2[i]);
      _off3 = ( (_off3<<1) + off3[i]);

      off = ((off&4088) | (_off0));
      off = ((off&4039) | ((_off1)<<3));
      off = ((off&3647) | ((_off2)<<6));
      off = ((off&511 ) | ((_off3)<<9));
    }
  }

  // CHECKING

  // the following trivial is_valid allows
  // the user of derived cell base classes
  // to add their own purpose checking

  bool is_valid(bool, int) const
  { return true; }

  // For use by Compact_container.
  void * for_compact_container() const { return N[0].for_compact_container(); }
  void for_compact_container(void *p) { N[0].for_compact_container(p); }

  // TDS internal data access functions.
  TDS_data& tds_data() { return _tds_data; }
  const TDS_data& tds_data() const { return _tds_data; }

  // TODO: Get rid of this flag! Used in convert_to_1_sheeted_covering.
  // Either use the conflict flag or a std::map.
  void set_additional_flag(unsigned char f) {
    CGAL_triangulation_assertion(f < 4);
    _additional_flag = f;
  }
  unsigned char get_additional_flag() const {
    return _additional_flag;
  }

private:
  Cell_handle   N[4];
  Vertex_handle V[4];
  TDS_data _tds_data;
  unsigned char _additional_flag:2;
  // 3 consecutive bits give the offset in x, y, and z.
  // From right to left:
  // - vertex(0): bit[0]-bit[2],
  // - vertex(1): bit[3]-bit[5],
  // - vertex(2): bit[6]-bit[8],
  // - vertex(3): bit[9]-bit[12].
  // Thus the underlying data type must have at least 12 bits,
  // which is true for uint.
  unsigned int off;
};

template < class TDS >
inline
std::istream&
operator>>(std::istream &is, Periodic_3_triangulation_ds_cell_base_3<TDS> &)
  // non combinatorial information. Default = nothing
{
  return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os,
    const Periodic_3_triangulation_ds_cell_base_3<TDS> &)
  // non combinatorial information. Default = nothing
{
  return os;
}

// Specialization for void.
template <>
class Periodic_3_triangulation_ds_cell_base_3<void>
{
public:
  typedef internal::Dummy_tds_3                  Triangulation_data_structure;
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Cell_handle     Cell_handle;
  template <typename TDS2>
  struct Rebind_TDS {
    typedef Periodic_3_triangulation_ds_cell_base_3<TDS2> Other;
  };
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_DS_CELL_BASE_3_H
