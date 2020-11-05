// Copyright (c) 2016  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

#ifndef CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_3_H
#define CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_3_H

#include <CGAL/license/Triangulation_3.h>


#include <CGAL/basic.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_ds_vertex_base_3.h>

namespace CGAL {

template < typename GT, typename DSVb = Triangulation_ds_vertex_base_3<> >
class Regular_triangulation_vertex_base_3
  : public DSVb
{
public:

  typedef typename DSVb::Cell_handle                   Cell_handle;
  typedef GT                                           Geom_traits;
  typedef typename GT::Weighted_point_3                Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename DSVb::template Rebind_TDS<TDS2>::Other  DSVb2;
    typedef Regular_triangulation_vertex_base_3<GT, DSVb2>           Other;
  };

  Regular_triangulation_vertex_base_3()
    : DSVb() {}

  Regular_triangulation_vertex_base_3(const Point & p)
    : DSVb(), _p(p) {}

  Regular_triangulation_vertex_base_3(const Point & p, Cell_handle c)
    : DSVb(c), _p(p) {}

  Regular_triangulation_vertex_base_3(Cell_handle c)
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
operator>>(std::istream &is, Regular_triangulation_vertex_base_3<GT, DSVb> &v)
  // non combinatorial information. Default = point
{
  return is >> static_cast<DSVb&>(v) >> v.point();
}

template < class GT, class DSVb >
std::ostream&
operator<<(std::ostream &os, const Regular_triangulation_vertex_base_3<GT, DSVb> &v)
  // non combinatorial information. Default = point
{
  return os << static_cast<const DSVb&>(v) << v.point();
}

} //namespace CGAL

#endif // CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_3_H
