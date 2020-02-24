// Copyright (c) 2001  
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
// Author(s)     : Sylvain Pion
//                 Mael Rouxel-Labbé

#ifndef CGAL_KERNEL_CHECKER_H
#define CGAL_KERNEL_CHECKER_H

// This file contains the definition of a kernel traits checker.

#include <CGAL/basic.h>
#include <CGAL/result_of.h>
#include <CGAL/use.h>

#include <utility>
#include <typeinfo>

namespace CGAL {

// Small utility to manipulate pairs for kernel objects, and
// simple things for bool, Sign...  Object is yet another case...
template < typename T1, typename T2 >
struct Pairify {
  typedef std::pair<T1, T2>  result_type;
  result_type operator()(const T1 &t1, const T2 &t2) const
  { return std::make_pair(t1, t2); }
};

template <>
struct Pairify <bool, bool> {
  typedef bool   result_type;
  result_type operator()(const bool &t1, const bool &t2) const
  { CGAL_kernel_assertion(t1 == t2); CGAL_USE(t2); return t1; }
};

template <>
struct Pairify <Sign, Sign> {
  typedef Sign   result_type;
  result_type operator()(const Sign &t1, const Sign &t2) const
  { CGAL_kernel_assertion(t1 == t2); CGAL_USE(t2); return t1; }
};

template <>
struct Pairify <Bounded_side, Bounded_side> {
  typedef Bounded_side   result_type;
  result_type operator()(const Bounded_side &t1, const Bounded_side &t2) const
  { CGAL_kernel_assertion(t1 == t2); CGAL_USE(t2); return t1; }
};

template <>
struct Pairify <Angle, Angle> {
  typedef Angle   result_type;
  result_type operator()(const Angle &t1, const Angle &t2) const
  { CGAL_kernel_assertion(t1 == t2); CGAL_USE(t2); return t1; }
};

// Class used by Kernel_checker.
template <class P1, class P2, class Cmp>
class Primitive_checker
{
  P1  p1;
  P2  p2;
  Cmp cmp;

public:
  Primitive_checker(const P1 &pp1 = P1(), const P2 &pp2 = P2(), const Cmp &c = Cmp())
    : p1(pp1), p2(pp2), cmp(c)
  { }

  template <class A1>
  typename Pairify<typename CGAL::cpp11::result_of<P1(const A1&)>::type,
                   typename CGAL::cpp11::result_of<P2(const A1&)>::type>::result_type
  operator()(const A1 &a1) const
  {
    typedef typename CGAL::cpp11::result_of<P1(const A1&)>::type result_type_1;
    typedef typename CGAL::cpp11::result_of<P2(const A1&)>::type result_type_2;
    result_type_1 res1 = p1(a1.first);
    result_type_2 res2 = p2(a1.second);
    if (! cmp(res1, res2))
    {
      std::cerr << "Kernel_checker error : " << res1 << " != " << res2
                << " for the inputs : " << std::endl;
      std::cerr << a1.first << std::endl;
      std::cerr << a1.second << std::endl;
      std::cerr << "functor first kernel : "
                << typeid(p1).name() << std::endl;
      std::cerr << "functor second kernel: "
                << typeid(p2).name() << std::endl;
      std::cerr << CGAL_PRETTY_FUNCTION << std::endl;
      CGAL_kernel_assertion(false);
    }
    return Pairify<result_type_1, result_type_2>()(res1, res2);
  }

  template <class A1, class A2>
  typename Pairify<typename CGAL::cpp11::result_of<P1(const A1&, const A2&)>::type,
                   typename CGAL::cpp11::result_of<P2(const A1&, const A2&)>::type>::result_type
  operator()(const A1 &a1, const A2 &a2) const
  {
    typedef typename CGAL::cpp11::result_of<P1(const A1&, const A2&)>::type result_type_1;
    typedef typename CGAL::cpp11::result_of<P2(const A1&, const A2&)>::type result_type_2;
    result_type_1 res1 = p1(a1.first, a2.first);
    result_type_2 res2 = p2(a1.second, a2.second);
    if (! cmp(res1, res2))
    {
      std::cerr << "Kernel_checker error : " << res1 << " != " << res2
                << " for the inputs : " << std::endl;
      std::cerr << a1.first << std::endl;
      std::cerr << a1.second << std::endl;
      std::cerr << a2.first << std::endl;
      std::cerr << a2.second << std::endl;
      std::cerr << "functor first kernel : "
                << typeid(p1).name() << std::endl;
      std::cerr << "functor second kernel: "
                << typeid(p2).name() << std::endl;
      std::cerr << CGAL_PRETTY_FUNCTION << std::endl;
      CGAL_kernel_assertion(false);
    }
    return Pairify<result_type_1, result_type_2>()(res1, res2);
  }

  template <class A1, class A2, class A3>
  typename Pairify<typename CGAL::cpp11::result_of<P1(
                     const A1&, const A2&, const A3&)>::type,
                   typename CGAL::cpp11::result_of<P2(
                     const A1&, const A2&, const A3&)>::type>::result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
  {
    typedef typename CGAL::cpp11::result_of<P1(const A1&, const A2&, const A3&)>::type result_type_1;
    typedef typename CGAL::cpp11::result_of<P2(const A1&, const A2&, const A3&)>::type result_type_2;
    result_type_1 res1 = p1(a1.first, a2.first, a3.first);
    result_type_2 res2 = p2(a1.second, a2.second, a3.second);
    if (! cmp(res1, res2))
    {
      std::cerr << "Kernel_checker error : " << res1 << " != " << res2
                << " for the inputs : " << std::endl;
      std::cerr << a1.first << std::endl;
      std::cerr << a1.second << std::endl;
      std::cerr << a2.first << std::endl;
      std::cerr << a2.second << std::endl;
      std::cerr << a3.first << std::endl;
      std::cerr << a3.second << std::endl;
      std::cerr << "functor first kernel : "
                << typeid(p1).name() << std::endl;
      std::cerr << "functor second kernel: "
                << typeid(p2).name() << std::endl;
      std::cerr << CGAL_PRETTY_FUNCTION << std::endl;
      CGAL_kernel_assertion(false);
    }
    return Pairify<result_type_1, result_type_2>()(res1, res2);
  }

  template <class A1, class A2, class A3, class A4>
  typename Pairify<typename CGAL::cpp11::result_of<P1(
                     const A1&, const A2&, const A3&, const A4&)>::type,
                   typename CGAL::cpp11::result_of<P2(
                     const A1&, const A2&, const A3&, const A4&)>::type>::result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
  {
    typedef typename CGAL::cpp11::result_of<P1(
                       const A1&, const A2&, const A3&, const A4&)>::type result_type_1;
    typedef typename CGAL::cpp11::result_of<P2(
                       const A1&, const A2&, const A3&, const A4&)>::type result_type_2;

    result_type_1 res1 = p1(a1.first, a2.first, a3.first, a4.first);
    result_type_2 res2 = p2(a1.second, a2.second, a3.second, a4.second);
    if (! cmp(res1, res2))
    {
      std::cerr << "Kernel_checker error : " << res1 << " != " << res2
                << " for the inputs : " << std::endl;
      std::cerr << a1.first << std::endl;
      std::cerr << a1.second << std::endl;
      std::cerr << a2.first << std::endl;
      std::cerr << a2.second << std::endl;
      std::cerr << a3.first << std::endl;
      std::cerr << a3.second << std::endl;
      std::cerr << a4.first << std::endl;
      std::cerr << a4.second << std::endl;
      std::cerr << "functor first kernel : "
                << typeid(p1).name() << std::endl;
      std::cerr << "functor second kernel: "
                << typeid(p2).name() << std::endl;
      std::cerr << CGAL_PRETTY_FUNCTION << std::endl;
      CGAL_kernel_assertion(false);
    }
    return Pairify<result_type_1, result_type_2>()(res1, res2);
  }

  template <class A1, class A2, class A3, class A4, class A5>
  typename Pairify<typename CGAL::cpp11::result_of<P1(
                     const A1&, const A2&, const A3&, const A4&, const A5&)>::type,
                   typename CGAL::cpp11::result_of<P2(
                     const A1&, const A2&, const A3&, const A4&, const A5&)>::type>::result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5) const
  {
    typedef typename CGAL::cpp11::result_of<P1(
      const A1&, const A2&, const A3&, const A4&, const A5&)>::type result_type_1;
    typedef typename CGAL::cpp11::result_of<P2(
      const A1&, const A2&, const A3&, const A4&, const A5&)>::type result_type_2;

    result_type_1 res1 = p1(a1.first, a2.first, a3.first, a4.first, a5.first);
    result_type_2 res2 = p2(a1.second, a2.second, a3.second, a4.second, a5.second);
    if (! cmp(res1, res2))
    {
      std::cerr << "Kernel_checker error : " << res1 << " != " << res2
                << " for the inputs : " << std::endl;
      std::cerr << a1.first << std::endl;
      std::cerr << a1.second << std::endl;
      std::cerr << a2.first << std::endl;
      std::cerr << a2.second << std::endl;
      std::cerr << a3.first << std::endl;
      std::cerr << a3.second << std::endl;
      std::cerr << a4.first << std::endl;
      std::cerr << a4.second << std::endl;
      std::cerr << a5.first << std::endl;
      std::cerr << a5.second << std::endl;
      std::cerr << "functor first kernel : "
                << typeid(p1).name() << std::endl;
      std::cerr << "functor second kernel: "
                << typeid(p2).name() << std::endl;
      std::cerr << CGAL_PRETTY_FUNCTION << std::endl;
      CGAL_kernel_assertion(false);
    }
    return Pairify<result_type_1, result_type_2>()(res1, res2);
  }

  // Same thing with more arguments...
};

struct dont_check_equal {
  template < typename T1, typename T2 >
  bool operator()(const T1 & /* t1 */, const T2 &/*t2*/) const
  { return true; }
  template < typename T >
  bool operator()(const T &t1, const T &t2) const
  { return t1 == t2; }
};

template < class K1, class K2, class Cmp = dont_check_equal >
class Kernel_checker
{
protected:
  K1 k1;
  K2 k2;
  Cmp cmp;

public:

  typedef bool                      Boolean;
  typedef CGAL::Sign                Sign;
  typedef CGAL::Comparison_result   Comparison_result;
  typedef CGAL::Orientation         Orientation;
  typedef CGAL::Oriented_side       Oriented_side;
  typedef CGAL::Bounded_side        Bounded_side;
  typedef CGAL::Angle               Angle;

  typedef K1     Kernel1;
  typedef K2     Kernel2;
  typedef Cmp    Comparator;

  // Kernel objects are defined as pairs, with primitives run in parallel.
#define CGAL_kc_pair(X) typedef std::pair<typename K1::X, typename K2::X> X;

  CGAL_kc_pair(RT)
  CGAL_kc_pair(FT)

  // TODO : Object_[23] are subtil : should probably be Object(pair<...>).
  // Or should Assign_[23] be used, and that's it ?
  // In any case, Assign will have to be treated separately because it
  // takes its first argument by non-const reference.
  // Maybe Primitive_checker should provide a variant with non-const ref...

  CGAL_kc_pair(Object_2)
  CGAL_kc_pair(Object_3)

  CGAL_kc_pair(Point_2)
  CGAL_kc_pair(Weighted_point_2)
  CGAL_kc_pair(Vector_2)
  CGAL_kc_pair(Direction_2)
  CGAL_kc_pair(Line_2)
  CGAL_kc_pair(Ray_2)
  CGAL_kc_pair(Segment_2)
  CGAL_kc_pair(Triangle_2)
  CGAL_kc_pair(Iso_rectangle_2)
  CGAL_kc_pair(Circle_2)
  CGAL_kc_pair(Conic_2)
  CGAL_kc_pair(Aff_transformation_2)

  CGAL_kc_pair(Point_3)
  CGAL_kc_pair(Weighted_point_3)
  CGAL_kc_pair(Plane_3)
  CGAL_kc_pair(Vector_3)
  CGAL_kc_pair(Direction_3)
  CGAL_kc_pair(Line_3)
  CGAL_kc_pair(Ray_3)
  CGAL_kc_pair(Segment_3)
  CGAL_kc_pair(Triangle_3)
  CGAL_kc_pair(Tetrahedron_3)
  CGAL_kc_pair(Iso_cuboid_3)
  CGAL_kc_pair(Sphere_3)
  CGAL_kc_pair(Aff_transformation_3)

#undef CGAL_kc_pair

#define CGAL_Kernel_pred(X, Y) \
  typedef Primitive_checker<typename K1::X, typename K2::X, Cmp> X; \
  X Y() const { return X(k1.Y(), k2.Y(), cmp); }

#define CGAL_Kernel_cons(Y,Z) CGAL_Kernel_pred(Y,Z)

  public:

  #include <CGAL/Kernel/interface_macros.h>
};

} //namespace CGAL

#endif // CGAL_KERNEL_CHECKER_H
