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

CGAL_BEGIN_NAMESPACE

// use to deduce arity of functors
// --> allows binding std functors
template < class T >
struct Arity_traits {
  enum { arity = T::arity };
};

// --------------------------------------------------------------------
// specializations for std functors:
//

template < class T >
struct Arity_traits< std::plus< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::minus< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::multiplies< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::divides< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::modulus< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::negate< T > > { enum { arity = 1 }; };
template < class T >
struct Arity_traits< std::equal_to< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::not_equal_to< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::greater< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::less< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::greater_equal< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::less_equal< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::logical_and< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::logical_or< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::logical_not< T > > { enum { arity = 1 }; };
template < class T >
struct Arity_traits< std::unary_negate< T > > { enum { arity = 1 }; };
template < class T >
struct Arity_traits< std::binary_negate< T > > { enum { arity = 2 }; };
template < class T >
struct Arity_traits< std::binder1st< T > > { enum { arity = 1 }; };
template < class T >
struct Arity_traits< std::binder2nd< T > > { enum { arity = 1 }; };
template < class T1, class T2 >
struct Arity_traits< std::pointer_to_unary_function< T1, T2 > > {
  enum { arity = 1 }; };
template < class T1, class T2, class T3 >
struct Arity_traits< std::pointer_to_binary_function< T1, T2, T3 > > {
  enum { arity = 2 }; };
template < class T1, class T2 >
struct Arity_traits< std::mem_fun_t< T1, T2 > > {
  enum { arity = 1 }; };
template < class T1, class T2, class T3 >
struct Arity_traits< std::mem_fun1_t< T1, T2, T3 > > {
  enum { arity = 2 }; };
template < class T1, class T2 >
struct Arity_traits< std::mem_fun_ref_t< T1, T2 > > {
  enum { arity = 1 }; };
template < class T1, class T2, class T3 >
struct Arity_traits< std::mem_fun1_ref_t< T1, T2, T3 > > {
  enum { arity = 2 }; };
template < class T1, class T2 >
struct Arity_traits< std::const_mem_fun_t< T1, T2 > > {
  enum { arity = 1 }; };
template < class T1, class T2, class T3 >
struct Arity_traits< std::const_mem_fun1_t< T1, T2, T3 > > {
  enum { arity = 2 }; };
template < class T1, class T2 >
struct Arity_traits< std::const_mem_fun_ref_t< T1, T2 > > {
  enum { arity = 1 }; };
template < class T1, class T2, class T3 >
struct Arity_traits< std::const_mem_fun1_ref_t< T1, T2, T3 > > {
  enum { arity = 2 }; };


template < class T >
struct Constant {
  enum { arity = 0 };
  Constant(const T& t_) : t(t_) {}
  T operator()() const { return t; }
protected:
  T t;
};
// Using (F::arity - 1) here gives an ICE on gcc 2.95
template < class F, int a, class A, int i >
struct Binder;

template < class F, class A >
struct Binder< F, 1, A, 1 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

  Binder(const F& f_, const A& a_) : f(f_), a(a_) {}

   
  result_type operator()( ) const
  { return f(a); }

protected:
  F f;
  A a;
};
template < class F, class A >
struct Binder< F, 2, A, 1 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

  Binder(const F& f_, const A& a_) : f(f_), a(a_) {}

  template < class A1 >
  result_type operator()(const A1& a1) const
  { return f(a, a1); }

protected:
  F f;
  A a;
};
template < class F, class A >
struct Binder< F, 2, A, 2 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

  Binder(const F& f_, const A& a_) : f(f_), a(a_) {}

  template < class A1 >
  result_type operator()(const A1& a1) const
  { return f(a1, a); }

protected:
  F f;
  A a;
};
template < class F, class A >
struct Binder< F, 3, A, 1 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

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
struct Binder< F, 3, A, 2 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

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
struct Binder< F, 3, A, 3 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

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
struct Binder< F, 4, A, 1 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

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
struct Binder< F, 4, A, 2 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

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
struct Binder< F, 4, A, 3 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

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
struct Binder< F, 4, A, 4 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

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
struct Binder< F, 5, A, 1 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

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
struct Binder< F, 5, A, 2 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

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
struct Binder< F, 5, A, 3 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

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
struct Binder< F, 5, A, 4 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

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
struct Binder< F, 5, A, 5 > {
  typedef typename F::result_type result_type;
  enum { arity = Arity_traits< F >::arity - 1 };

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

template < class T, class A, int i >
struct Binder_class {
  typedef Binder< T, Arity_traits< T >::arity, A, i > Type;
};
template < class F, class A >
inline typename Binder_class< F, A, 1 >::Type
bind_1(const F& f, const A& a) {
  typedef typename Binder_class< F, A, 1 >::Type B;
  return B(f, a);
}

template < class F, class A >
inline typename Binder_class< F, A, 2 >::Type
bind_2(const F& f, const A& a) {
  typedef typename Binder_class< F, A, 2 >::Type B;
  return B(f, a);
}

template < class F, class A >
inline typename Binder_class< F, A, 3 >::Type
bind_3(const F& f, const A& a) {
  typedef typename Binder_class< F, A, 3 >::Type B;
  return B(f, a);
}

template < class F, class A >
inline typename Binder_class< F, A, 4 >::Type
bind_4(const F& f, const A& a) {
  typedef typename Binder_class< F, A, 4 >::Type B;
  return B(f, a);
}

template < class F, class A >
inline typename Binder_class< F, A, 5 >::Type
bind_5(const F& f, const A& a) {
  typedef typename Binder_class< F, A, 5 >::Type B;
  return B(f, a);
}


namespace CGALi {
  struct Not_used { enum { arity = -1 }; };
}

template < class F0,
           int a0,
           class F1,
           int a1 = Arity_traits< F1 >::arity,
           class F2 = CGALi::Not_used,
           int a2 = Arity_traits< F2 >::arity,
           class F3 = CGALi::Not_used,
           int a3 = Arity_traits< F3 >::arity >
struct Compose;

// ------------------------------------------------------------------------
// one function to compose
// ------------------------------------------------------------------------

template < class F0, class F1 >
struct Compose< F0, 1,
                F1, 0,
                CGALi::Not_used, -1,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 0 };

  Compose(const F0& f0_, const F1& f1_) : f0(f0_), f1(f1_) {}

   
  result_type operator()( ) const
  { return f0(f1( )); }

protected:
  F0 f0;
  F1 f1;
};
template < class F0, class F1 >
struct Compose< F0, 1,
                F1, 1,
                CGALi::Not_used, -1,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 1 };

  Compose(const F0& f0_, const F1& f1_) : f0(f0_), f1(f1_) {}

  template < class A1 >
  result_type operator()(const A1& a1) const
  { return f0(f1(a1)); }

protected:
  F0 f0;
  F1 f1;
};
template < class F0, class F1 >
struct Compose< F0, 1,
                F1, 2,
                CGALi::Not_used, -1,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 2 };

  Compose(const F0& f0_, const F1& f1_) : f0(f0_), f1(f1_) {}

  template < class A1, class A2 >
  result_type operator()(const A1& a1, const A2& a2) const
  { return f0(f1(a1, a2)); }

protected:
  F0 f0;
  F1 f1;
};
template < class F0, class F1 >
struct Compose< F0, 1,
                F1, 3,
                CGALi::Not_used, -1,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 3 };

  Compose(const F0& f0_, const F1& f1_) : f0(f0_), f1(f1_) {}

  template < class A1, class A2, class A3 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3) const
  { return f0(f1(a1, a2, a3)); }

protected:
  F0 f0;
  F1 f1;
};
template < class F0, class F1 >
struct Compose< F0, 1,
                F1, 4,
                CGALi::Not_used, -1,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 4 };

  Compose(const F0& f0_, const F1& f1_) : f0(f0_), f1(f1_) {}

  template < class A1, class A2, class A3, class A4 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4) const
  { return f0(f1(a1, a2, a3, a4)); }

protected:
  F0 f0;
  F1 f1;
};
template < class F0, class F1 >
struct Compose< F0, 1,
                F1, 5,
                CGALi::Not_used, -1,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 5 };

  Compose(const F0& f0_, const F1& f1_) : f0(f0_), f1(f1_) {}

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

// unary functions
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 1,
                F2, 0,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 1 + 0 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 0,
                F2, 1,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 0 + 1 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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

// binary functions
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 2,
                F2, 0,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 2 + 0 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 1,
                F2, 1,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 1 + 1 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 0,
                F2, 2,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 0 + 2 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 3,
                F2, 0,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 3 + 0 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 2,
                F2, 1,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 2 + 1 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 1,
                F2, 2,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 1 + 2 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 0,
                F2, 3,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 0 + 3 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 4,
                F2, 0,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 4 + 0 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 3,
                F2, 1,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 3 + 1 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 2,
                F2, 2,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 2 + 2 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 1,
                F2, 3,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 1 + 3 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 0,
                F2, 4,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 0 + 4 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 5,
                F2, 0,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 5 + 0 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 4,
                F2, 1,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 4 + 1 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 3,
                F2, 2,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 3 + 2 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 2,
                F2, 3,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 2 + 3 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 1,
                F2, 4,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 1 + 4 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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
struct Compose< F0, 2,
                F1, 0,
                F2, 5,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 0 + 5 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
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

// 6-arg functions
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 5,
                F2, 1,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 5 + 1 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6) const
  { return f0(f1(a1, a2, a3, a4, a5), f2(a6)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 4,
                F2, 2,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 4 + 2 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6) const
  { return f0(f1(a1, a2, a3, a4), f2(a5, a6)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 3,
                F2, 3,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 3 + 3 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6) const
  { return f0(f1(a1, a2, a3), f2(a4, a5, a6)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 2,
                F2, 4,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 2 + 4 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6) const
  { return f0(f1(a1, a2), f2(a3, a4, a5, a6)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 1,
                F2, 5,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 1 + 5 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6) const
  { return f0(f1(a1), f2(a2, a3, a4, a5, a6)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};

// 7-arg functions
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 5,
                F2, 2,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 5 + 2 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6, class A7 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6, const A7& a7) const
  { return f0(f1(a1, a2, a3, a4, a5), f2(a6, a7)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 4,
                F2, 3,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 4 + 3 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6, class A7 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6, const A7& a7) const
  { return f0(f1(a1, a2, a3, a4), f2(a5, a6, a7)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 3,
                F2, 4,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 3 + 4 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6, class A7 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6, const A7& a7) const
  { return f0(f1(a1, a2, a3), f2(a4, a5, a6, a7)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 2,
                F2, 5,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 2 + 5 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6, class A7 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6, const A7& a7) const
  { return f0(f1(a1, a2), f2(a3, a4, a5, a6, a7)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};

// 8-arg functions
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 5,
                F2, 3,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 5 + 3 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6, class A7, class A8 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6, const A7& a7, const A8& a8) const
  { return f0(f1(a1, a2, a3, a4, a5), f2(a6, a7, a8)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 4,
                F2, 4,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 4 + 4 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6, class A7, class A8 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6, const A7& a7, const A8& a8) const
  { return f0(f1(a1, a2, a3, a4), f2(a5, a6, a7, a8)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 3,
                F2, 5,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 3 + 5 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6, class A7, class A8 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6, const A7& a7, const A8& a8) const
  { return f0(f1(a1, a2, a3), f2(a4, a5, a6, a7, a8)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};

// 9-arg functions
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 5,
                F2, 4,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 5 + 4 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6, class A7, class A8, class A9 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6, const A7& a7, const A8& a8,
                         const A9& a9) const
  { return f0(f1(a1, a2, a3, a4, a5), f2(a6, a7, a8, a9)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 4,
                F2, 5,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 4 + 5 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6, class A7, class A8, class A9 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6, const A7& a7, const A8& a8,
                         const A9& a9) const
  { return f0(f1(a1, a2, a3, a4), f2(a5, a6, a7, a8, a9)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};

// 10-arg functions
template < class F0, class F1, class F2 >
struct Compose< F0, 2,
                F1, 5,
                F2, 5,
                CGALi::Not_used, -1 >
{
  typedef typename F0::result_type result_type;
  enum { arity = 5 + 5 };

  Compose(const F0& f0_, const F1& f1_, const F2& f2_)
  : f0(f0_), f1(f1_), f2(f2_)
  {}

  template < class A1, class A2, class A3, class A4, class A5,
             class A6, class A7, class A8, class A9, class A10 >
  result_type operator()(const A1& a1, const A2& a2,
                         const A3& a3, const A4& a4, const A5& a5,
                         const A6& a6, const A7& a7, const A8& a8,
                         const A9& a9, const A10& a10) const
  { return f0(f1(a1, a2, a3, a4), f2(a5, a6, a7, a8, a9, a10)); }

protected:
  F0 f0;
  F1 f1;
  F2 f2;
};


template < class F0, class F1,
           class F2 = CGALi::Not_used,
           class F3 = CGALi::Not_used >
struct Compose_class {
  typedef Compose< F0, Arity_traits< F0 >::arity,
                   F1, Arity_traits< F1 >::arity,
                   F2, Arity_traits< F2 >::arity > Type;
};

template < class F0, class F1 >
inline typename Compose_class< F0, F1 >::Type
compose(const F0& f0, const F1& f1) {
  typedef typename Compose_class< F0, F1 >::Type C;
  return C(f0, f1);
}

template < class F0, class F1, class F2 >
inline typename Compose_class< F0, F1, F2 >::Type
compose(const F0& f0, const F1& f1, const F2& f2)
{
  typedef typename Compose_class< F0, F1, F2 >::Type C;
  return C(f0, f1, f2);
}

template < class F0, class F1, class F2, class F3 >
inline typename Compose_class< F0, F1, F2, F3 >::Type
compose(const F0& f0, const F1& f1, const F2& f2, const F3& f3)
{
  typedef typename Compose_class< F0, F1, F2, F3 >::Type C;
  return C(f0, f1, f2, f3);
}



CGAL_END_NAMESPACE

#endif // CGAL_FUNCTIONAL_H //
// EOF //
