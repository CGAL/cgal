// ============================================================================
//
// Copyright (c) 1997, 1998, 1999, 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : functional.h
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// package       : $CGAL_Package: STL_Extension $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@cs.unc.edu>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// New Functor Adaptors.
// ============================================================================

#ifndef CGAL_FUNCTIONAL_H
#define CGAL_FUNCTIONAL_H 1

#include <CGAL/functional_base.h>

#ifdef CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION
#include <CGAL/functional_msvc.h>
#else

CGAL_BEGIN_NAMESPACE

// +----------------------------------------------------------------------+
// | Functor Adaptors for Swapping Arguments
// +----------------------------------------------------------------------+

namespace CGALi {

  template < class F, int i, class A > struct Swapper;

  template < class F >
  struct Swapper< F, 1, Arity_tag< 2 > > {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper(const F& f_) : f(f_) {}
  
    template < class A1, class A2 >
    result_type operator()
    (const A1& a1, const A2& a2) const
    { return f(a2, a1); }
  
  protected:
    F f;
  };
  template < class F >
  struct Swapper< F, 1, Arity_tag< 3 > > {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper(const F& f_) : f(f_) {}
  
    template < class A1, class A2, class A3 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3) const
    { return f(a2, a1, a3); }
  
  protected:
    F f;
  };
  template < class F >
  struct Swapper< F, 2, Arity_tag< 3 > > {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper(const F& f_) : f(f_) {}
  
    template < class A1, class A2, class A3 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3) const
    { return f(a1, a3, a2); }
  
  protected:
    F f;
  };
  template < class F >
  struct Swapper< F, 1, Arity_tag< 4 > > {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper(const F& f_) : f(f_) {}
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return f(a2, a1, a3, a4); }
  
  protected:
    F f;
  };
  template < class F >
  struct Swapper< F, 2, Arity_tag< 4 > > {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper(const F& f_) : f(f_) {}
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return f(a1, a3, a2, a4); }
  
  protected:
    F f;
  };
  template < class F >
  struct Swapper< F, 3, Arity_tag< 4 > > {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper(const F& f_) : f(f_) {}
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return f(a1, a2, a4, a3); }
  
  protected:
    F f;
  };
  template < class F >
  struct Swapper< F, 1, Arity_tag< 5 > > {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper(const F& f_) : f(f_) {}
  
    template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5) const
    { return f(a2, a1, a3, a4, a5); }
  
  protected:
    F f;
  };
  template < class F >
  struct Swapper< F, 2, Arity_tag< 5 > > {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper(const F& f_) : f(f_) {}
  
    template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5) const
    { return f(a1, a3, a2, a4, a5); }
  
  protected:
    F f;
  };
  template < class F >
  struct Swapper< F, 3, Arity_tag< 5 > > {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper(const F& f_) : f(f_) {}
  
    template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5) const
    { return f(a1, a2, a4, a3, a5); }
  
  protected:
    F f;
  };
  template < class F >
  struct Swapper< F, 4, Arity_tag< 5 > > {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper(const F& f_) : f(f_) {}
  
    template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5) const
    { return f(a1, a2, a3, a5, a4); }
  
  protected:
    F f;
  };

} // namespace CGALi

template < class F, int i = 1 >
struct Swap {
  typedef CGALi::Swapper< F, i, typename Arity_traits< F >::Arity > Type;
};

template < class F > inline
typename Swap< F, 1 >::Type
swap_1(const F& f) {
  typedef typename Swap< F, 1 >::Type S;
  return S(f);
}

template < class F > inline
typename Swap< F, 2 >::Type
swap_2(const F& f) {
  typedef typename Swap< F, 2 >::Type S;
  return S(f);
}

template < class F > inline
typename Swap< F, 3 >::Type
swap_3(const F& f) {
  typedef typename Swap< F, 3 >::Type S;
  return S(f);
}

template < class F > inline
typename Swap< F, 4 >::Type
swap_4(const F& f) {
  typedef typename Swap< F, 4 >::Type S;
  return S(f);
}

// +----------------------------------------------------------------------+
// | Binding Arguments of a Functor
// +----------------------------------------------------------------------+

namespace CGALi {
  // Using (F::arity - 1) here gives an ICE on gcc 2.95
  template < class F, class R, class A, int i >
  struct Binder;

  template < class F, class A >
  struct Binder< F, Arity_tag< 1 >, A, 1 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 1 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
     
    result_type operator()( ) const
    { return f(a); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 2 >, A, 1 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 2 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1 >
    result_type operator()(const A1& a1) const
    { return f(a, a1); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 2 >, A, 2 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 2 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1 >
    result_type operator()(const A1& a1) const
    { return f(a1, a); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 3 >, A, 1 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 3 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1, class A2 >
    result_type operator()(
                             const A1& a1, const A2& a2) const
    { return f(a, a1, a2); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 3 >, A, 2 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 3 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1, class A2 >
    result_type operator()(
                             const A1& a1, const A2& a2) const
    { return f(a1, a, a2); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 3 >, A, 3 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 3 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1, class A2 >
    result_type operator()(
                             const A1& a1, const A2& a2) const
    { return f(a1, a2, a); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 4 >, A, 1 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 4 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f(a, a1, a2, a3); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 4 >, A, 2 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 4 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f(a1, a, a2, a3); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 4 >, A, 3 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 4 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f(a1, a2, a, a3); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 4 >, A, 4 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 4 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f(a1, a2, a3, a); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 5 >, A, 1 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 5 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f(a, a1, a2, a3, a4); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 5 >, A, 2 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 5 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f(a1, a, a2, a3, a4); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 5 >, A, 3 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 5 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f(a1, a2, a, a3, a4); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 5 >, A, 4 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 5 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f(a1, a2, a3, a, a4); }
  
  protected:
    F f;
    A a;
  };
  template < class F, class A >
  struct Binder< F, Arity_tag< 5 >, A, 5 > {
    typedef typename F::result_type  result_type;
    typedef Arity_tag< 5 - 1 >      Arity;
  
    Binder(const F& f_, const A& a_) : f(f_), a(a_) {}
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f(a1, a2, a3, a4, a); }
  
  protected:
    F f;
    A a;
  };

} // namespace CGALi

template < class T, class A, int i >
struct Bind {
  typedef CGALi::Binder< T, typename Arity_traits< T >::Arity, A, i > Type;
};
template < class F, class A >
inline typename Bind< F, A, 1 >::Type
bind_1(const F& f, const A& a) {
  typedef typename Bind< F, A, 1 >::Type B;
  return B(f, a);
}

template < class F, class A >
inline typename Bind< F, A, 2 >::Type
bind_2(const F& f, const A& a) {
  typedef typename Bind< F, A, 2 >::Type B;
  return B(f, a);
}

template < class F, class A >
inline typename Bind< F, A, 3 >::Type
bind_3(const F& f, const A& a) {
  typedef typename Bind< F, A, 3 >::Type B;
  return B(f, a);
}

template < class F, class A >
inline typename Bind< F, A, 4 >::Type
bind_4(const F& f, const A& a) {
  typedef typename Bind< F, A, 4 >::Type B;
  return B(f, a);
}

template < class F, class A >
inline typename Bind< F, A, 5 >::Type
bind_5(const F& f, const A& a) {
  typedef typename Bind< F, A, 5 >::Type B;
  return B(f, a);
}


// +----------------------------------------------------------------------+
// | Composing Functors
// +----------------------------------------------------------------------+

namespace CGALi {
  struct Not_used { typedef Arity_tag< -1 > Arity; };

  template < class F0, class A0,
             class F1, class A1,
             class F2 = CGALi::Not_used, class A2 = Arity_tag< -1 >,
             class F3 = CGALi::Not_used, class A3 = Arity_tag< -1 > >
  struct Composer;

  template < class F0, class A0,
             class F1, class A1,
             class F2, class A2,
             class F3 = CGALi::Not_used, class A3 = Arity_tag< -1 > >
  struct Composer_shared;

  // ------------------------------------------------------------------------
  // one function to compose
  // ------------------------------------------------------------------------

  template < class F0, class F1 >
  struct Composer< F0, Arity_tag< 1 >,
                   F1, Arity_tag< 0 >,
                   CGALi::Not_used, Arity_tag< -1 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_) : f0(f0_), f1(f1_) {}
  
     
    result_type operator()( ) const
    { return f0(f1( )); }
  
  protected:
    F0 f0;
    F1 f1;
  };
  template < class F0, class F1 >
  struct Composer< F0, Arity_tag< 1 >,
                   F1, Arity_tag< 1 >,
                   CGALi::Not_used, Arity_tag< -1 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_) : f0(f0_), f1(f1_) {}
  
    template < class A1 >
    result_type operator()(const A1& a1) const
    { return f0(f1(a1)); }
  
  protected:
    F0 f0;
    F1 f1;
  };
  template < class F0, class F1 >
  struct Composer< F0, Arity_tag< 1 >,
                   F1, Arity_tag< 2 >,
                   CGALi::Not_used, Arity_tag< -1 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_) : f0(f0_), f1(f1_) {}
  
    template < class A1, class A2 >
    result_type operator()(const A1& a1, const A2& a2) const
    { return f0(f1(a1, a2)); }
  
  protected:
    F0 f0;
    F1 f1;
  };
  template < class F0, class F1 >
  struct Composer< F0, Arity_tag< 1 >,
                   F1, Arity_tag< 3 >,
                   CGALi::Not_used, Arity_tag< -1 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 3 > Arity;
  
    Composer(const F0& f0_, const F1& f1_) : f0(f0_), f1(f1_) {}
  
    template < class A1, class A2, class A3 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3) const
    { return f0(f1(a1, a2, a3)); }
  
  protected:
    F0 f0;
    F1 f1;
  };
  template < class F0, class F1 >
  struct Composer< F0, Arity_tag< 1 >,
                   F1, Arity_tag< 4 >,
                   CGALi::Not_used, Arity_tag< -1 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 4 > Arity;
  
    Composer(const F0& f0_, const F1& f1_) : f0(f0_), f1(f1_) {}
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()(
                             const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4) const
    { return f0(f1(a1, a2, a3, a4)); }
  
  protected:
    F0 f0;
    F1 f1;
  };
  template < class F0, class F1 >
  struct Composer< F0, Arity_tag< 1 >,
                   F1, Arity_tag< 5 >,
                   CGALi::Not_used, Arity_tag< -1 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 5 > Arity;
  
    Composer(const F0& f0_, const F1& f1_) : f0(f0_), f1(f1_) {}
  
    template < class A1, class A2, class A3, class A4,
      class A5 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2, a3, a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
  };

  // ------------------------------------------------------------------------
  // two functions to compose
  // ------------------------------------------------------------------------

  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 0 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
     
    result_type operator()( ) const
    { return f0(f1( ), f2( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };

  // unary functions
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 0 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    template < class A1 >
    result_type operator()(const A1& a1) const
    { return f0(f1(a1), f2( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 1 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    template < class A1 >
    result_type operator()(const A1& a1) const
    { return f0(f1( ), f2(a1)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };

  // binary functions
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 0 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    template < class A1, class A2 >
    result_type operator()(const A1& a1, const A2& a2) const
    { return f0(f1(a1, a2), f2( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 1 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    template < class A1, class A2 >
    result_type operator()(const A1& a1, const A2& a2) const
    { return f0(f1(a1), f2(a2)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 2 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    template < class A1, class A2 >
    result_type operator()(const A1& a1, const A2& a2) const
    { return f0(f1(), f2(a1, a2)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };

  // 3-arg functions
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 3 >,
                   F2, Arity_tag< 0 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 3 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    template < class A1, class A2, class A3 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3) const
    { return f0(f1(a1, a2, a3), f2( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 1 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    template < class A1, class A2, class A3 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3) const
    { return f0(f1(a1, a2), f2(a3)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 2 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    template < class A1, class A2, class A3 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3) const
    { return f0(f1(a1), f2(a2, a3)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 3 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 3 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    template < class A1, class A2, class A3 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3) const
    { return f0(f1( ), f2(a1, a2, a3)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };

  // 4-arg functions
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 4 >,
                   F2, Arity_tag< 0 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 4 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    
      template < class A1, class A2, class A3, class A4 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4) const
    { return f0(f1(a1, a2, a3, a4), f2( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 3 >,
                   F2, Arity_tag< 1 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 3 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    
      template < class A1, class A2, class A3, class A4 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4) const
    { return f0(f1(a1, a2, a3), f2(a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 2 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    
      template < class A1, class A2, class A3, class A4 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4) const
    { return f0(f1(a1, a2), f2(a3, a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 3 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 3 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    
      template < class A1, class A2, class A3, class A4 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4) const
    { return f0(f1(a1), f2(a2, a3, a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 4 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 4 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    
      template < class A1, class A2, class A3, class A4 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4) const
    { return f0(f1( ), f2(a1, a2, a3, a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };

  // 5-arg functions
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 5 >,
                   F2, Arity_tag< 0 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 5 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    
      template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2, a3, a4, a5), f2( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 4 >,
                   F2, Arity_tag< 1 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 4 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    
      template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2, a3, a4), f2(a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 3 >,
                   F2, Arity_tag< 2 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 3 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    
      template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2, a3), f2(a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 3 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 3 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    
      template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2), f2(a3, a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 4 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 4 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    
      template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4, const A5& a5) const
    { return f0(f1(a1), f2(a2, a3, a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer< F0, Arity_tag< 2 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 5 >,
                   CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 5 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    
      template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,
                             const A3& a3, const A4& a4, const A5& a5) const
    { return f0(f1( ), f2(a1, a2, a3, a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };

                                          
                                  
                          
                      
            
  template < class F0, class F1, class F2 >
  struct Composer_shared< F0, Arity_tag< 2 >,
                          F1, Arity_tag< 0 >,
                          F2, Arity_tag< 0 >,
                          CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 > Arity;
  
    Composer_shared(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
     
    result_type
    operator()( ) const
    { return f0(f1( ), f2( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer_shared< F0, Arity_tag< 2 >,
                          F1, Arity_tag< 1 >,
                          F2, Arity_tag< 1 >,
                          CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 > Arity;
  
    Composer_shared(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    template < class A1 >
    result_type
    operator()(const A1& a1) const
    { return f0(f1(a1), f2(a1)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer_shared< F0, Arity_tag< 2 >,
                          F1, Arity_tag< 2 >,
                          F2, Arity_tag< 2 >,
                          CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 > Arity;
  
    Composer_shared(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
    template < class A1, class A2 >
    result_type
    operator()(const A1& a1, const A2& a2) const
    { return f0(f1(a1, a2), f2(a1, a2)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer_shared< F0, Arity_tag< 2 >,
                          F1, Arity_tag< 3 >,
                          F2, Arity_tag< 3 >,
                          CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 3 > Arity;
  
    Composer_shared(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
      template < class A1, class A2, class A3 >
    result_type
    operator()(  const A1& a1, const A2& a2, const A3& a3) const
    { return f0(f1(  a1, a2, a3), f2(  a1, a2, a3)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer_shared< F0, Arity_tag< 2 >,
                          F1, Arity_tag< 4 >,
                          F2, Arity_tag< 4 >,
                          CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 4 > Arity;
  
    Composer_shared(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
      template < class A1, class A2, class A3, class A4 >
    result_type
    operator()(  const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return f0(f1(  a1, a2, a3, a4), f2(  a1, a2, a3, a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };
  template < class F0, class F1, class F2 >
  struct Composer_shared< F0, Arity_tag< 2 >,
                          F1, Arity_tag< 5 >,
                          F2, Arity_tag< 5 >,
                          CGALi::Not_used, Arity_tag< -1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 5 > Arity;
  
    Composer_shared(const F0& f0_, const F1& f1_, const F2& f2_)
    : f0(f0_), f1(f1_), f2(f2_)
    {}
  
      template < class A1, class A2, class A3, class A4, class A5 >
    result_type
    operator()(  const A1& a1, const A2& a2, const A3& a3,
                 const A4& a4, const A5& a5) const
    { return f0(f1(  a1, a2, a3, a4, a5), f2(  a1, a2, a3, a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
  };


  // ------------------------------------------------------------------------
  // three functions to compose
  // ------------------------------------------------------------------------

  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 0 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
     
    result_type operator()( ) const
    { return f0(f1( ), f2( ), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };

  // unary functions
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 0 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    template < class A1 >
    result_type operator()(const A1& a1) const
    { return f0(f1(a1), f2( ), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 1 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    template < class A1 >
    result_type operator()(const A1& a1) const
    { return f0(f1( ), f2(a1), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 0 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    template < class A1 >
    result_type operator()(const A1& a1) const
    { return f0(f1( ), f2( ), f3(a1)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };

  // binary functions
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 0 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2 >
    result_type operator()(const A1& a1, const A2& a2) const
    { return f0(f1(a1, a2), f2( ), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 2 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 2 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2 >
    result_type operator()(const A1& a1, const A2& a2) const
    { return f0(f1( ), f2(a1, a2), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 2 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 0 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2 >
    result_type operator()(const A1& a1, const A2& a2) const
    { return f0(f1( ), f2( ), f3(a1, a2)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 1 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2 >
    result_type operator()(const A1& a1, const A2& a2) const
    { return f0(f1(a1), f2(a2), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 0 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2 >
    result_type operator()(const A1& a1, const A2& a2) const
    { return f0(f1(a1), f2( ), f3(a2 )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 1 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2 >
    result_type operator()(const A1& a1, const A2& a2) const
    { return f0(f1( ), f2(a1), f3(a2 )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };

  // 3-arg functions
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 3 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 3 + 0 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f0(f1(a1, a2, a3), f2( ), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 3 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 3 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f0(f1( ), f2(a1, a2, a3), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 3 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 0 + 3 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f0(f1( ), f2( ), f3(a1, a2, a3)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 1 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f0(f1(a1, a2), f2(a3), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 0 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f0(f1(a1, a2), f2( ), f3(a3)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 2 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 2 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f0(f1(a1), f2(a2, a3), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 2 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 2 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f0(f1( ), f2(a1, a2), f3(a3)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 2 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 0 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f0(f1(a1), f2( ), f3(a2, a3)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 2 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 1 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f0(f1( ), f2(a1), f3(a2, a3)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 1 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3 >
    result_type operator()(
                             const A1& a1, const A2& a2, const A3& a3) const
    { return f0(f1(a1), f2(a2), f3(a3)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };

  // 4-arg functions
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 4 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 4 + 0 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1(a1, a2, a3, a4), f2( ), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 4 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 4 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1( ), f2(a1, a2, a3, a4), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 4 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 0 + 4 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1( ), f2( ), f3(a1, a2, a3, a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 3 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 3 + 1 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1(a1, a2, a3), f2(a4), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 3 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 3 + 0 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1(a1, a2, a3), f2( ), f3(a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 3 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 3 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1(a1), f2(a2, a3, a4), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 3 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 3 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1( ), f2(a1, a2, a3), f3(a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 3 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 0 + 3 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1(a1), f2( ), f3(a2, a3, a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 3 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 1 + 3 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1( ), f2(a1), f3(a2, a3, a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 2 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 2 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1(a1, a2), f2(a3, a4), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 2 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 0 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1(a1, a2), f2( ), f3(a3, a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 2 >,
                   F3, Arity_tag< 2 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 2 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1( ), f2(a1, a2), f3(a3, a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 1 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1(a1, a2), f2(a3), f3(a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 2 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 2 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1(a1), f2(a2, a3), f3(a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 2 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 1 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4 >
    result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                             const A4& a4) const
    { return f0(f1(a1), f2(a2), f3(a3, a4)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };

  // 5-arg functions
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 5 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 5 + 0 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2, a3, a4, a5), f2( ), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 5 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 5 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1( ), f2(a1, a2, a3, a4, a5), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 5 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 0 + 5 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1( ), f2( ), f3(a1, a2, a3, a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 4 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 4 + 1 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2, a3, a4), f2(a5), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 4 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 4 + 0 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2, a3, a4), f2( ), f3(a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 4 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 4 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1), f2(a2, a3, a4, a5), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 4 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 4 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1( ), f2(a1, a2, a3, a4), f3(a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 4 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 0 + 4 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1), f2( ), f3(a2, a3, a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 4 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 1 + 4 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1( ), f2(a1), f3(a2, a3, a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 3 >,
                   F2, Arity_tag< 2 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 3 + 2 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2, a3), f2(a4, a5), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 3 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 2 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 3 + 0 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2, a3), f2( ), f3(a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 3 >,
                   F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 3 + 0 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2), f2(a3, a4, a5), f3( )); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 3 >,
                   F3, Arity_tag< 2 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 3 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1( ), f2(a1, a2, a3), f3(a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 0 >,
                   F3, Arity_tag< 3 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 0 + 3 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2), f2( ), f3(a3, a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 0 >,
                   F2, Arity_tag< 2 >,
                   F3, Arity_tag< 3 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 + 2 + 3 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1( ), f2(a1, a2), f3(a3, a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 3 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 3 + 1 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2, a3), f2(a4), f3(a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 3 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 3 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1), f2(a2, a3, a4), f3(a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 3 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 1 + 3 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1), f2(a2), f3(a3, a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 2 >,
                   F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 2 + 1 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2), f2(a3, a4), f3(a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 2 >,
                   F2, Arity_tag< 1 >,
                   F3, Arity_tag< 2 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 + 1 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1, a2), f2(a3), f3(a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer< F0, Arity_tag< 3 >,
                   F1, Arity_tag< 1 >,
                   F2, Arity_tag< 2 >,
                   F3, Arity_tag< 2 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 + 2 + 2 > Arity;
  
    Composer(const F0& f0_, const F1& f1_, const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    
      template < class A1, class A2, class A3,
      class A4, class A5 >
    result_type operator()(const A1& a1, const A2& a2,const A3& a3,
                             const A4& a4, const A5& a5) const
    { return f0(f1(a1), f2(a2, a3), f3(a4, a5)); }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };

  template < class F0, class F1, class F2, class F3 >
  struct Composer_shared< F0, Arity_tag< 3 >,
                          F1, Arity_tag< 0 >,
                          F2, Arity_tag< 0 >,
                          F3, Arity_tag< 0 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 0 > Arity;
  
    Composer_shared(const F0& f0_, const F1& f1_,
                    const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
     
    result_type
    operator()( ) const
    { return f0(f1( ),
                f2( ),
                f3( ));
    }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer_shared< F0, Arity_tag< 3 >,
                          F1, Arity_tag< 1 >,
                          F2, Arity_tag< 1 >,
                          F3, Arity_tag< 1 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 1 > Arity;
  
    Composer_shared(const F0& f0_, const F1& f1_,
                    const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    template < class A1 >
    result_type
    operator()(const A1& a1) const
    { return f0(f1(a1),
                f2(a1),
                f3(a1));
    }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer_shared< F0, Arity_tag< 3 >,
                          F1, Arity_tag< 2 >,
                          F2, Arity_tag< 2 >,
                          F3, Arity_tag< 2 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 2 > Arity;
  
    Composer_shared(const F0& f0_, const F1& f1_,
                    const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
    template < class A1, class A2 >
    result_type
    operator()(const A1& a1, const A2& a2) const
    { return f0(f1(a1, a2),
                f2(a1, a2),
                f3(a1, a2));
    }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer_shared< F0, Arity_tag< 3 >,
                          F1, Arity_tag< 3 >,
                          F2, Arity_tag< 3 >,
                          F3, Arity_tag< 3 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 3 > Arity;
  
    Composer_shared(const F0& f0_, const F1& f1_,
                    const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
      template < class A1, class A2, class A3 >
    result_type
    operator()(  const A1& a1, const A2& a2, const A3& a3) const
    { return f0(f1(  a1, a2, a3),
                f2(  a1, a2, a3),
                f3(  a1, a2, a3));
    }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer_shared< F0, Arity_tag< 3 >,
                          F1, Arity_tag< 4 >,
                          F2, Arity_tag< 4 >,
                          F3, Arity_tag< 4 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 4 > Arity;
  
    Composer_shared(const F0& f0_, const F1& f1_,
                    const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
      template < class A1, class A2, class A3, class A4 >
    result_type
    operator()(  const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return f0(f1(  a1, a2, a3, a4),
                f2(  a1, a2, a3, a4),
                f3(  a1, a2, a3, a4));
    }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };
  template < class F0, class F1, class F2, class F3 >
  struct Composer_shared< F0, Arity_tag< 3 >,
                          F1, Arity_tag< 5 >,
                          F2, Arity_tag< 5 >,
                          F3, Arity_tag< 5 > >
  {
    typedef typename F0::result_type result_type;
    typedef Arity_tag< 5 > Arity;
  
    Composer_shared(const F0& f0_, const F1& f1_,
                    const F2& f2_, const F3& f3_)
    : f0(f0_), f1(f1_), f2(f2_), f3(f3_)
    {}
  
      template < class A1, class A2, class A3, class A4, class A5 >
    result_type
    operator()(  const A1& a1, const A2& a2, const A3& a3,
                 const A4& a4, const A5& a5) const
    { return f0(f1(  a1, a2, a3, a4, a5),
                f2(  a1, a2, a3, a4, a5),
                f3(  a1, a2, a3, a4, a5));
    }
  
  protected:
    F0 f0;
    F1 f1;
    F2 f2;
    F3 f3;
  };

} // namespace CGALi

// ------------------------------------------------------------------------
// Encapsulate the composition type ==> can be adapted
// ------------------------------------------------------------------------

template < class F0, class F1,
           class F2 = CGALi::Not_used,
           class F3 = CGALi::Not_used >
struct Compose {
  typedef CGALi::Composer< F0, typename Arity_traits< F0 >::Arity,
                           F1, typename Arity_traits< F1 >::Arity,
                           F2, typename Arity_traits< F2 >::Arity,
                           F3, typename Arity_traits< F3 >::Arity > Type;
};

template < class F0, class F1, class F2, class F3 = CGALi::Not_used >
struct Compose_shared {
  typedef CGALi::Composer_shared< F0, typename Arity_traits< F0 >::Arity,
                                  F1, typename Arity_traits< F1 >::Arity,
                                  F2, typename Arity_traits< F2 >::Arity,
                                  F3, typename Arity_traits< F3 >::Arity >
  Type;
};
// ------------------------------------------------------------------------
// compose helper functions
// ------------------------------------------------------------------------

template < class F0, class F1 >
inline typename Compose< F0, F1 >::Type
compose(const F0& f0, const F1& f1) {
  typedef typename Compose< F0, F1 >::Type C;
  return C(f0, f1);
}

template < class F0, class F1, class F2 >
inline typename Compose< F0, F1, F2 >::Type
compose(const F0& f0, const F1& f1, const F2& f2)
{
  typedef typename Compose< F0, F1, F2 >::Type C;
  return C(f0, f1, f2);
}

template < class F0, class F1, class F2 >
inline typename Compose_shared< F0, F1, F2 >::Type
compose_shared(const F0& f0, const F1& f1, const F2& f2)
{
  typedef typename Compose_shared< F0, F1, F2 >::Type C;
  return C(f0, f1, f2);
}

template < class F0, class F1, class F2, class F3 >
inline typename Compose< F0, F1, F2, F3 >::Type
compose(const F0& f0, const F1& f1, const F2& f2, const F3& f3)
{
  typedef typename Compose< F0, F1, F2, F3 >::Type C;
  return C(f0, f1, f2, f3);
}

template < class F0, class F1, class F2, class F3 >
inline typename Compose_shared< F0, F1, F2, F3 >::Type
compose_shared(const F0& f0, const F1& f1, const F2& f2, const F3& f3)
{
  typedef typename Compose_shared< F0, F1, F2, F3 >::Type C;
  return C(f0, f1, f2, f3);
}

CGAL_END_NAMESPACE

#endif // CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION

#endif // CGAL_FUNCTIONAL_H //
// EOF //
