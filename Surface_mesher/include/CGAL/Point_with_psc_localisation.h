// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POINT_WITH_PSC_LOCALISATION_H
#define CGAL_POINT_WITH_PSC_LOCALISATION_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Point_traits.h>

#include <string>

namespace CGAL {

template <class Point, typename Index_type = int>
class Point_with_psc_localisation : public Point
{
  typedef CGAL::Point_traits<Point> Point_traits;
  typedef typename Point_traits::Bare_point Bare_point;
  typedef typename Kernel_traits<Bare_point>::Kernel Kernel;
  typedef typename Kernel::FT FT;
  typedef Point_with_psc_localisation<Point> Self;

public:
  typedef Tag_true Set_on_surface_tag;

  Point_with_psc_localisation() : Point(), index(), dim(-1) {}

  Point_with_psc_localisation(const Point& p) : Point(p), index(), dim(-1) {}

//   Point_with_psc_localisation(const typename Kernel::Point_3& p) : Point(p), index(0) {}

  Point_with_psc_localisation(const Point_with_psc_localisation& p)
    : Point(p), index(p.element_index()), dim(p.dimension()) {}

  Point_with_psc_localisation(const FT& x, const FT& y, const FT& z, const FT& w = FT(1))
    : Point(Point_traits().point(Bare_point(x, y, z, w))), index(), dim(-1) {}

  const Index_type& element_index() const
  {
    return index;
  }

  void set_element_index(const Index_type i)
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

  void set_on_surface(const Index_type i)
  {
//     CGAL_assertion(i>=0);
    set_element_index(i);
    set_dimension(2);
  }

  void set_on_curve(const Index_type i)
  {
//     CGAL_assertion(i>=0);
    set_element_index(i);
    set_dimension(1);
  }

  void set_on_vertex(const Index_type i)
  {
//     CGAL_assertion(i>=0);
    set_element_index(i);
    set_dimension(0);
  }

#ifdef CGAL_MESH_3_IO_H
  static
  std::string io_signature()
  {
    return Get_io_signature<Point>()() + "+i+i";
  }
#endif
private:
  Index_type index;
  int dim;
}; // end class Point_with_psc_localisation

template <class Point>
class Point_traits<Point_with_psc_localisation<Point> >
  : public Point_traits<Point>
{
};

template <class Point>
struct Is_weighted<Point_with_psc_localisation<Point> >
  : public Is_weighted<Point>
{
};

template <class Point>
std::ostream&
operator<<(std::ostream &os, const Point_with_psc_localisation<Point>& p)
{
  os << static_cast<const Point&>(p);
  if(is_ascii(os))
    os << ' ' << p.dimension() << ' ' << p.element_index();
  else {
    write(os, p.dimension());
    write(os, p.element_index());
  }
  return os;
}

template <class Point>
std::istream&
operator>>(std::istream &is, Point_with_psc_localisation<Point>& p)
{
  is >>  static_cast<Point&>(p);
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
