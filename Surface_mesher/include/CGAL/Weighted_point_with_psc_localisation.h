// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
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
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POINT_WITH_PSC_LOCALISATION_H
#define CGAL_POINT_WITH_PSC_LOCALISATION_H

#include <CGAL/Point_traits.h>

#include <string>

namespace CGAL {

template <class Weighted_point>
class Weighted_point_with_psc_localisation : public Weighted_point
{
  typedef CGAL::Point_traits<Weighted_point> Point_traits;
  typedef typename Point_traits::Bare_point Bare_point;
  typedef typename Kernel_traits<Bare_point>::Kernel Kernel;
  typedef typename Kernel::FT FT;

public:
  Weighted_point_with_psc_localisation() : Weighted_point(), index(-1), dim(-1) {}

  Weighted_point_with_psc_localisation(const Weighted_point& p) : Weighted_point(p), index(-1), dim(-1) {}
			   
  Weighted_point_with_psc_localisation(const Bare_point& p) : Weighted_point(p), index(-1), dim(-1) {}
			   
//   Weighted_point_with_psc_localisation(const typename Kernel::Point_3& p) : Weighted_point(p), index(-1), dim(-1) {}

  Weighted_point_with_psc_localisation(const typename Kernel::Point_3& p,
                                       const FT& w) : Weighted_point(p, w), index(-1), dim(-1) {}

  Weighted_point_with_psc_localisation(const Weighted_point_with_psc_localisation& p)
    : Weighted_point(p), index(p.element_index()), dim(p.dimension()) {}

  Weighted_point_with_psc_localisation(const FT& x, const FT& y, const FT& z, const FT& w = FT(1))
    : Weighted_point(Point_traits().point(Bare_point(x, y, z, w))), index(-1), dim(-1) {}

  const int& element_index() const
  {
    return index;
  }

  void set_element_index(const int i)
  {
    index = i;
  }

  void set_dimension(const int d)
  {
    dim = d;
  }

  const int& dimension() const
  {
    return dim;
  }

  void set_on_surface(const int i)
  {
    CGAL_assertion(i>=0);
    set_element_index(i);
    set_dimension(2);
  }

  void set_on_curve(const int i)
  {
    CGAL_assertion(i>=0);
    set_element_index(i);
    set_dimension(1);
  }

  void set_on_vertex(const int i)
  {
    CGAL_assertion(i>=0);
    set_element_index(i);
    set_dimension(0);
  }

#ifdef CGAL_MESH_3_IO_H
  static
  std::string io_signature()
  {
    return Get_io_signature<Weighted_point>()() + "+i+i";
  }
#endif
private:
  int index;
  int dim;
}; // end class Weighted_point_with_psc_localisation

template <class Weighted_point>
class Point_traits<Weighted_point_with_psc_localisation<Weighted_point> >
  : public Point_traits<Weighted_point>
{
};

template <class Weighted_point>
class Is_weighted<Weighted_point_with_psc_localisation<Weighted_point> >
  : public Is_weighted<Weighted_point>
{
};

template <class Weighted_point>
std::ostream&
operator<<(std::ostream &os, const Weighted_point_with_psc_localisation<Weighted_point>& p)
{
  os << static_cast<const Weighted_point&>(p);
  if(is_ascii(os))
    os << ' ' << p.dimension() << ' ' << p.element_index();
  else {
    write(os, p.dimension());
    write(os, p.element_index());
  }
  return os;
}

template <class Weighted_point>
std::istream&
operator>>(std::istream &is, Weighted_point_with_psc_localisation<Weighted_point>& p)
{
  is >>  static_cast<Weighted_point&>(p);
  int index, dim;
  if(is_ascii(is))
    is >> dim >> index;
  else {
    read(is, dim);
    read(is, index);
  }
  p.set_dimension(dim);
  p.set_element_index(index);
  return is;
}

} // end namespace CGAL

#endif
