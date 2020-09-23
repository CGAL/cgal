// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_KERNEL_D_FUNCTION_OBJECTS_H
#define CGAL_KERNEL_D_FUNCTION_OBJECTS_H

#include <CGAL/intersections_d.h>

// These functors come from the 2D-3D kernels.
// Since they have changed there, they now need to be copied here.

namespace CGAL {

namespace internal {

template <class ToBeConstructed>
class Construct
{
  public:
    typedef ToBeConstructed  result_type;

    template <class ... A>
    ToBeConstructed
    operator()( A&& ... a) const
    { return ToBeConstructed(std::forward<A>(a)...); }
};

class Call_has_on_positive_side
{
  public:
    typedef bool           result_type;

    template <class Cls, class Arg>
    bool
    operator()( const Cls& c, const Arg& a) const
    { return c.has_on_positive_side(a); }
};

class Call_oriented_side
{
  public:
    typedef Oriented_side   result_type;

    template <class Cls, class Arg>
    Oriented_side
    operator()( const Cls& c, const Arg& a) const
    { return c.oriented_side(a); }
};

template<class R>
class Intersect
{
private:
  typedef typename R::Line_d       Line_d;
  typedef typename R::Point_d      Point_d;
  typedef typename R::Segment_d    Segment_d;
  typedef typename R::Ray_d        Ray_d;
  typedef typename R::Hyperplane_d Hyperplane_d;

public:
  template <typename>
  struct result;

  template <typename F>
  struct result<F(Line_d, Line_d)>
  { typedef boost::optional< boost::variant< Point_d, Line_d > > type; };

  template <typename F>
  struct result<F(Segment_d, Line_d)>
  { typedef boost::optional< boost::variant< Point_d, Segment_d > > type; };
  template <typename F>
  struct result<F(Line_d, Segment_d)> : result<F(Segment_d, Line_d)>
  { };

  template <typename F>
  struct result<F(Segment_d, Segment_d)>
  { typedef boost::optional< boost::variant< Point_d, Segment_d > > type; };

  template <typename F>
  struct result<F(Ray_d, Line_d)>
  { typedef boost::optional< boost::variant< Point_d, Ray_d > > type; };

  template <typename F>
  struct result<F(Line_d, Ray_d)> : result<F(Ray_d, Line_d)>
  { };

  template <typename F>
  struct result<F(Ray_d, Segment_d)>
  { typedef boost::optional< boost::variant< Point_d, Segment_d > > type; };

  template <typename F>
  struct result<F(Segment_d, Ray_d)> : result<F(Ray_d, Segment_d)>
  { };

  template <typename F>
  struct result<F(Ray_d, Ray_d)>
  { typedef boost::optional< boost::variant< Point_d, Segment_d, Ray_d > > type; };

  template <typename F>
  struct result<F(Hyperplane_d, Line_d)>
  { typedef boost::optional< boost::variant< Point_d, Line_d > > type; };
  template <typename F>
  struct result<F(Line_d, Hyperplane_d)> : result<F(Hyperplane_d, Line_d)>
  { };

  template <typename F>
  struct result<F(Hyperplane_d, Ray_d)>
  { typedef boost::optional< boost::variant< Point_d, Ray_d > > type; };
  template <typename F>
  struct result<F(Ray_d, Hyperplane_d)> : result<F(Hyperplane_d, Ray_d)>
  { };

  template <typename F>
  struct result<F(Hyperplane_d, Segment_d)>
  { typedef boost::optional< boost::variant< Point_d, Segment_d > > type; };
  template <typename F>
  struct result<F(Segment_d, Hyperplane_d)> : result<F(Hyperplane_d, Segment_d)>
  { };

  template <class T1, class T2>
  typename result<Intersect(T1,T2)>::type
  operator()(const T1& t1, const T2& t2) const
  { return Intersections::internal::intersection(t1, t2, R()); }
};

template<class R>
class Do_intersect
{
  public:
    typedef bool result_type;

    template <class T1, class T2>
    bool
    operator()(const T1& t1, const T2& t2) const
  { return CGAL::Intersections::internal::do_intersect(t1, t2, R()); }
};

} // end namespace internal
} //namespace CGAL

#endif // CGAL_KERNEL_D_FUNCTION_OBJECTS_H
