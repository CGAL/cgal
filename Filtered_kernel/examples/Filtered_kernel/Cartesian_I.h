// Copyright (c) 2000,2001,2002,2003
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
// Author(s)     : Herve Bronnimann, Sylvain Pion

#ifndef CGAL_CARTESIAN_I_H
#define CGAL_CARTESIAN_I_H

#include <CGAL/Cartesian/Cartesian_base.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>

namespace CGAL {

namespace CartesianKernelFunctors {

  template <typename K>
  class Intersect_with_iterators_2
  {
public:
    typedef typename K::Point_2       Point_2;
    typedef typename K::Segment_2     Segment_2;
  public:
    typedef void result_type;

    Intersect_with_iterators_2() {}

    template <typename OutputIterator>
    OutputIterator
    operator()(const Segment_2& s1, const Segment_2& s2, OutputIterator it ) const
    {
      Object o = typename K::Intersect_2()(s1, s2);
      if(! o.is_empty()){
        *it = o;
        it++;
      }

      return it;
    }
  };

} // namespace CartesianKernelFunctors


template < typename FT_, typename Kernel_ >
struct Cartesian_base_ref_count_I
  : public Cartesian_base< Kernel_, FT_ >
{
    typedef FT_                                           RT;
    typedef FT_                                           FT;

    // The mechanism that allows to specify reference-counting or not.
    template < typename T >
    struct Handle { typedef Handle_for<T>    type; };

    template < typename Kernel2 >
    struct Base { typedef Cartesian_base_ref_count_I<FT_, Kernel2>  Type; };

    // TODO: cleanup (use Rational_traits<> instead)
    static   FT make_FT(const RT & num, const RT& denom) { return num/denom;}
    static   FT make_FT(const RT & num)                  { return num;}
    static   RT FT_numerator(const FT &r)                { return r;}
    static   RT FT_denominator(const FT &)               { return RT(1);}

  // was in Cartesian_base
  typedef Kernel_ K;


#define CGAL_Kernel_pred(Y,Z) typedef CartesianKernelFunctors::Y<K> Y; \
                              Y Z() const { return Y(); }
#define CGAL_Kernel_cons(Y,Z) CGAL_Kernel_pred(Y,Z)


  CGAL_Kernel_cons(Intersect_with_iterators_2,
                   intersect_with_iterators_2_object)
#include <CGAL/Kernel/interface_macros.h>

    template<typename A, typename B>
    struct Intersection_result : CGAL::Intersection_traits<K, A, B> {};
};

template < typename FT_ >
struct Cartesian_I
  : public Type_equality_wrapper<
                Cartesian_base_ref_count_I<FT_, Cartesian_I<FT_> >,
                Cartesian_I<FT_> >
{};

} //namespace CGAL

#endif // CGAL_CARTESIAN_I_H
