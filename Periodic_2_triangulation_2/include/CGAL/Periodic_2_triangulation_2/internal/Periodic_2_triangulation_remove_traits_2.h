// Copyright (c) 2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_2_TRIANGULATION_REMOVE_TRAITS_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_REMOVE_TRAITS_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Periodic_2_offset_2.h>

namespace CGAL {

// Triangulation_2 uses Construct_point_2 to handle weighted and bare points.
// The default Construct_point_2 inherited by Periodic_2_triangulation_remove_traits_2
// must be overwritten by a custom Construct_point_2 that offers:
// - pair<K::Point_2, offset> --> pair<K::Point_2, offset> (identity)
template<class Gt_, typename Construct_point_2_base_>
class Construct_point_from_pair_2
  : public Construct_point_2_base_
{
  typedef Construct_point_2_base_        Base;
  typedef Gt_                            Geom_traits;

  // Gt::Point_2 is actually a pair <K::Point_2, offset>
  typedef typename Geom_traits::Point_2  Point_2;

public:
  Construct_point_from_pair_2(const Base& cp) : Base(cp) { }

  using Base::operator();

  template<typename F>
  struct result : Base::template result<F> {};

  template<typename F>
  struct result<F(Point_2)> {
    typedef const Point_2& type;
  };

  const Point_2& operator()(const Point_2& p) const { return p; }
};

template < class Traits_, class Functor_ >
class Functor_with_point_offset_pair_adaptor
  : public Functor_
{
  typedef Traits_                        Traits;
  typedef Functor_                       Functor;

  // `Traits::Point_2` is actually a `std::pair<Point_2, Offset>`
  typedef typename Traits::Point_2       Point;

public:
  typedef typename Functor::result_type result_type;

  Functor_with_point_offset_pair_adaptor(const Functor & functor) : Functor_(functor) { }

public:
  using Functor::operator();

  result_type operator()(const Point& p0, const Point& p1) const
  {
    return operator()(p0.first, p1.first,
                      p0.second, p1.second);
  }
  result_type operator()(const Point& p0, const Point& p1, const Point& p2) const
  {
    return operator()(p0.first, p1.first, p2.first,
                      p0.second, p1.second, p2.second);
  }
  result_type operator()(const Point& p0, const Point& p1, const Point& p2, const Point& p3) const
  {
    return operator()(p0.first, p1.first, p2.first, p3.first,
                      p0.second, p1.second, p2.second, p3.second);
  }
};

template <class Gt_,
          class Off_ = typename CGAL::Periodic_2_offset_2>
class Periodic_2_triangulation_remove_traits_2
  : public Gt_
{
  typedef Periodic_2_triangulation_remove_traits_2<Gt_, Off_>           Self;
  typedef Gt_                                                           Base;

public:
  typedef Gt_                                                           Geom_traits;

  typedef Off_                                                          Offset;

  typedef typename Geom_traits::RT                                      RT;
  typedef typename Geom_traits::FT                                      FT;
  typedef typename Geom_traits::Point_2                                 Bare_point;
  typedef std::pair<Bare_point, Offset>                                 Point_2;
  typedef typename Geom_traits::Domain                                  Domain;

  Periodic_2_triangulation_remove_traits_2(const Geom_traits& gt) : Base(gt) { }

  // Construct point
  typedef Construct_point_from_pair_2<Self, typename Geom_traits::Construct_point_2> Construct_point_2;

  // Triangulation predicates
  typedef Functor_with_point_offset_pair_adaptor<Self, typename Geom_traits::Compare_xy_2>
      Compare_xy_2;
  typedef Functor_with_point_offset_pair_adaptor<Self, typename Geom_traits::Orientation_2>
      Orientation_2;

  // Operations
  Construct_point_2 construct_point_2_object() const {
    return Construct_point_2();
  }
  Compare_xy_2 compare_xy_2_object() const {
    return Compare_xy_2(this->Base::compare_xy_2_object());
  }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_REMOVE_TRAITS_2_H
