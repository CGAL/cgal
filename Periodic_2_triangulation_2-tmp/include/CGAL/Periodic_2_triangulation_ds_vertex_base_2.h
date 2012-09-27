// Copyright (c) 1999-2003,2007-2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://nicokruithof@scm.gforge.inria.fr/svnroot/cgal/trunk/Periodic_2_triangulation_2/include/CGAL/Periodic_2_triangulation_ds_vertex_base_2.h $
// $Id: Periodic_2_triangulation_ds_vertex_base_2.h 56667 2010-06-09 07:37:13Z sloriot $
// 
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_TRIANGULATION_DS_VERTEX_BASE_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_DS_VERTEX_BASE_2_H

#include <CGAL/basic.h>
#include <CGAL/Dummy_tds_2.h>
#include <CGAL/Periodic_2_offset_2.h>

namespace CGAL {


/* A vertex class with an additionnal handle */
template < class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
class Periodic_2_triangulation_ds_vertex_base_2
  : public  Vb
{
  typedef Vb                              Base;
public:
  typedef typename Vb::Vertex_handle      Vertex_handle;
  typedef typename Vb::Face_handle        Face_handle;
  typedef typename Vb::Point              Point;
  typedef Periodic_2_offset_2       Offset;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other     Vb2;
    typedef Periodic_2_triangulation_ds_vertex_base_2<Gt,Vb2> Other;
  };

public:
  Periodic_2_triangulation_ds_vertex_base_2() : Base() {}
  Periodic_2_triangulation_ds_vertex_base_2(const Point & p)
    : Base(p), _off(), _offset_flag(false) {}
  Periodic_2_triangulation_ds_vertex_base_2(const Point & p, Face_handle f)
    : Base(f,p), _off(), _offset_flag(false) {}
  Periodic_2_triangulation_ds_vertex_base_2(Face_handle f)
    : Base(f), _off(), _offset_flag(false) {}

  const Offset& offset() const {
    return _off;
  }
  
  void set_offset(const Offset& off) {
    _off = off; _offset_flag=true;
  }

  void clear_offset() {
    _offset_flag=false;
    _off = Offset();
  }

  bool get_offset_flag() const {
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


  // // For use by the Compact_container.
  // void *   for_compact_container() const
  // { return _c.for_compact_container(); }
  // void * & for_compact_container()
  // { return _c.for_compact_container(); }

// private:
//   Face_handle _c;
//   Offset _off;
//   int _index;
//   bool _offset_flag;
};

template < class TDS >
inline
std::istream&
operator>>(std::istream &is, Periodic_2_triangulation_ds_vertex_base_2<TDS> &)
  // no combinatorial information.
{
  return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os,
    const Periodic_2_triangulation_ds_vertex_base_2<TDS> &)
  // no combinatorial information.
{
  return os;
}

// Specialization for void.
template <>
class Periodic_2_triangulation_ds_vertex_base_2<void>
{
public:
  typedef Dummy_tds_2                                   Triangulation_data_structure;
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Face_handle     Face_handle;
  template <typename TDS2>
  struct Rebind_TDS {
    typedef Periodic_2_triangulation_ds_vertex_base_2<TDS2> Other;
  };
};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_DS_VERTEX_BASE_2_H
