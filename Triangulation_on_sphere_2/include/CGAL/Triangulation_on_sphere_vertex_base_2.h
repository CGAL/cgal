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
// Author(s)     : Mariette Yvinec

#ifndef CGAL_TRIANGULATION_ON_SPHERE_VERTEX_BASE_2_H
#define CGAL_TRIANGULATION_ON_SPHERE_VERTEX_BASE_2_H

#include <CGAL/license/Triangulation_on_sphere_2.h>

#include <CGAL/Triangulation_ds_vertex_base_2.h>

namespace CGAL {

template <typename GT,
          typename Vb = Triangulation_ds_vertex_base_2<> >
class Triangulation_on_sphere_vertex_base_2
  : public Vb
{
  typedef typename Vb::Triangulation_data_structure                  Tds;
public:
  typedef GT                                                         Geom_traits;
  typedef typename GT::Point_on_sphere_2                             Point;
  typedef Tds                                                        Triangulation_data_structure;
  typedef typename Tds::Face_handle                                  Face_handle;
  typedef typename Tds::Vertex_handle                                Vertex_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other            Vb2;
    typedef Triangulation_on_sphere_vertex_base_2<GT, Vb2>           Other;
  };

private:
  Point _p;

public:
  Triangulation_on_sphere_vertex_base_2 () : Vb(), _p() { }
  Triangulation_on_sphere_vertex_base_2(const Point & p) : Vb(), _p(p) { }
  Triangulation_on_sphere_vertex_base_2(const Point & p, Face_handle f) : Vb(f), _p(p) { }
  Triangulation_on_sphere_vertex_base_2(Face_handle f) : Vb(f) { }

  void set_point(const Point & p) { _p = p; }
  const Point& point() const { return _p; }

  // the non-const version of point() is undocument but needed to make the point iterator works
  // using Lutz projection scheme
  Point& point() { return _p; }

  // the following trivial is_valid to allow the user of derived face base classes
  // to add their own purpose checking
  bool is_valid(bool /* verbose */ = false, int /* level */ = 0) const {return true;}
};

template < class GT, class Vb >
std::istream&
operator>>(std::istream &is, Triangulation_on_sphere_vertex_base_2<GT, Vb> &v)
{
  // non combinatorial information. Default = point
  return is >> static_cast<Vb&>(v) >> v.point();
}

template < class GT, class Vb >
std::ostream&
operator<<(std::ostream &os, const Triangulation_on_sphere_vertex_base_2<GT, Vb> &v)
{
  // non combinatorial information. Default = point
  return os << static_cast<const Vb&>(v) << v.point();
}

} // namespace CGAL

#endif //CGAL_TRIANGULATION_ON_SPHERE_VERTEX_BASE_2_H
