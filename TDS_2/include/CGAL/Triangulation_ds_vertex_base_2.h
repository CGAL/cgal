// Copyright (c) 1999,2000,2001,2002,2003  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Mariette Yvinec


#ifndef CGAL_TRIANGULATION_DS_VERTEX_BASE_2_H
#define CGAL_TRIANGULATION_DS_VERTEX_BASE_2_H

#include <CGAL/license/TDS_2.h>


#include <CGAL/config.h>
#include <iostream>
#include <CGAL/Dummy_tds_2.h>

namespace CGAL {

template < class TDS = void >
class Triangulation_ds_vertex_base_2 
{

public:
  typedef TDS                          Triangulation_data_structure;
  typedef typename TDS::Face_handle    Face_handle;
  typedef typename TDS::Vertex_handle  Vertex_handle;

  template <typename TDS2>
  struct Rebind_TDS { typedef Triangulation_ds_vertex_base_2<TDS2> Other; };

  Triangulation_ds_vertex_base_2 ()    : _f()  {}
  Triangulation_ds_vertex_base_2(Face_handle f)    :  _f(f)    {}

  Face_handle face() const { return _f;}
  void set_face(Face_handle f) { _f = f ;}

  //the following trivial is_valid to allow
  // the user of derived face base classes 
  // to add their own purpose checking
  bool is_valid(bool /*verbose*/=false, int /*level*/= 0) const
    {return face() != Face_handle();}

    // For use by the Compact_container.
  void *   for_compact_container() const { return _f.for_compact_container(); }
  void * & for_compact_container()       { return _f.for_compact_container(); }

private:
  Face_handle _f;
};

// Specialization for void.
template <>
class Triangulation_ds_vertex_base_2<void>
{
public:
  typedef Dummy_tds_2  Triangulation_data_structure;
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Face_handle     Face_handle;
  template <typename TDS2>
  struct Rebind_TDS { typedef Triangulation_ds_vertex_base_2<TDS2> Other; };
};


template < class TDS >
inline
std::istream&
operator>>(std::istream &is, Triangulation_ds_vertex_base_2<TDS> &)
  // no combinatorial information.
{
  return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os, const Triangulation_ds_vertex_base_2<TDS> &)
  // no combinatorial information.
{
  return os;
}
} //namespace CGAL

#endif //CGAL_TRIANGULATION_DS_VERTEX_BASE_2_H
