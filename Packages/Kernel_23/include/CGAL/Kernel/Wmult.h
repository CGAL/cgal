// ============================================================================
//
// Copyright (c) 1998,2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Kernel/Wmult.h
// author(s)     : Geert-Jan Giezeman, Sylvain Pion
//
// coordinator   : Saarbruecken
//
// ============================================================================

#ifndef CGAL_KERNEL_WMULT_H
#define CGAL_KERNEL_WMULT_H

CGAL_BEGIN_NAMESPACE

template < typename Rep_Tag >
struct wmult_tag;

template <>
struct wmult_tag<Cartesian_tag>
{
  template < typename RT >
  const RT & operator()(const RT &a, const RT &) const
  { return a; }

  template < typename RT >
  const RT & operator()(const RT &a, const RT &, const RT &) const
  { return a; }

  template < typename RT >
  const RT & operator()(const RT &a, const RT &, const RT &, const RT &) const
  { return a; }

  template < typename RT >
  const RT & operator()(const RT &a, const RT &, const RT &, const RT &,
                        const RT &) const
  { return a; }
};

template <>
struct wmult_tag<Homogeneous_tag>
{
  template < typename RT >
  RT operator()(const RT &a, const RT &w) const
  { return a*w; }

  template < typename RT >
  RT operator()(const RT &a, const RT &w1, const RT &w2) const
  { return a*w1*w2; }

  template < typename RT >
  RT operator()(const RT &a, const RT &w1, const RT &w2, const RT &w3) const
  { return a*w1*w2*w3; }

  template < typename RT >
  RT operator()(const RT &a, const RT &w1, const RT &w2, const RT &w3,
                const RT &w4) const
  { return a*w1*w2*w3*w4; }
};

template < typename K >
struct wmult_functor
  : wmult_tag<typename K::Rep_tag> {};

template < typename K, typename RT >
inline
RT wmult(K*, const RT &a, const RT &w)
{
    return wmult_functor<K>()(a, w);
}

template < typename K, typename RT >
inline
RT wmult(K*, const RT &a, const RT &w1, const RT &w2)
{
    return wmult_functor<K>()(a, w1, w2);
}

template < typename K, typename RT >
inline
RT wmult(K*, const RT &a, const RT &w1, const RT &w2, const RT &w3)
{
    return wmult_functor<K>()(a, w1, w2, w3);
}

template < typename K, typename RT >
inline
RT
wmult(K*, const RT &a, const RT &w1, const RT &w2, const RT &w3, const RT &w4)
{
    return wmult_functor<K>()(a, w1, w2, w3, w4);
}

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_WMULT_H
