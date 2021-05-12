// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri


#ifndef CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_2_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/Triangulation_ds_vertex_base_2.h>

namespace CGAL {

template < class GT,
           class Vb = Triangulation_ds_vertex_base_2<> >
class Regular_triangulation_vertex_base_2
  :   public Vb
{
  typedef typename Vb::Triangulation_data_structure     TDS;

public:
  typedef GT                                  Geom_traits;
  typedef typename GT::Point_2                Bare_point;
  typedef typename GT::Weighted_point_2       Weighted_point;
  typedef Weighted_point                      Point;
  typedef TDS                                 Triangulation_data_structure;
  typedef typename TDS::Face_handle           Face_handle;
  typedef typename TDS::Vertex_handle         Vertex_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other     Vb2;
    typedef Regular_triangulation_vertex_base_2<GT, Vb2>      Other;
  };

private:
  Weighted_point _p;
  bool _hidden;

public:
  Regular_triangulation_vertex_base_2 ()
    : Vb(), _hidden(false) {}
  Regular_triangulation_vertex_base_2(const Weighted_point & p)
    :  Vb(), _p(p), _hidden(false) {}
  Regular_triangulation_vertex_base_2(const Weighted_point & p, Face_handle f)
    : Vb(f),_p(p), _hidden(false) {}

  void set_hidden(bool b) { _hidden = b; }
  bool is_hidden() { return _hidden ;}
  void set_point(const Weighted_point & p) { _p = p; }
  const Weighted_point&  point() const { return _p; }
  Weighted_point&  point() { return _p; }
  bool is_valid(bool /* verbose */ = false, int /* level */ = 0) const
    {return true;}
};



template < class GT, class Vb >
std::istream&
operator>>(std::istream &is,
           Regular_triangulation_vertex_base_2<GT, Vb> &v)
  // non combinatorial information. Default = point
{
  return is >> static_cast<Vb&>(v) >> v.point();
}

template < class GT, class Vb >
std::ostream&
operator<<(std::ostream &os,
           const Regular_triangulation_vertex_base_2<GT, Vb> &v)
  // non combinatorial information. Default = point
{
  return os << static_cast<const Vb&>(v) << v.point();
}

} //namespace CGAL


#endif //CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_2_H
