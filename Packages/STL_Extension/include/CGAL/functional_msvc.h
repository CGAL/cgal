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
// file          : functional_msvc.h
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
// New Functor Adaptors (not using partial spec.).
// ============================================================================

#ifndef CGAL_FUNCTIONAL_MSVC_H
#define CGAL_FUNCTIONAL_MSVC_H 1

CGAL_BEGIN_NAMESPACE

template < class F, class A >
struct Binder_1 {
  typedef typename F::result_type result_type;
  Binder_1(const F& f_, const A& a_) : f(f_), a(a_) {}

  result_type
  operator()() const
  { return f(a); }

  template < class A1 >
  result_type
  operator()(const A1& a1) const
  { return f(a, a1); }

  template < class A1, class A2 >
  result_type
  operator()(const A1& a1, const A2& a2) const
  { return f(a, a1, a2); }

  template < class A1, class A2, class A3 >
  result_type
  operator()(const A1& a1, const A2& a2, const A3& a3) const
  { return f(a, a1, a2, a3); }

  template < class A1, class A2, class A3, class A4 >
  result_type
  operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
  { return f(a, a1, a2, a3, a4); }

protected:
  F f;
  A a;
};

template < class F, class A >
struct Binder_2 {
  typedef typename F::result_type result_type;
  Binder_2(const F& f_, const A& a_) : f(f_), a(a_) {}

  template < class A1 >
  result_type
  operator()(const A1& a1) const
  { return f(a1, a); }

  template < class A1, class A2 >
  result_type
  operator()(const A1& a1, const A2& a2) const
  { return f(a1, a, a2); }

  template < class A1, class A2, class A3 >
  result_type
  operator()(const A1& a1, const A2& a2, const A3& a3) const
  { return f(a1, a, a2, a3); }

  template < class A1, class A2, class A3, class A4 >
  result_type
  operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
  { return f(a1, a, a2, a3, a4); }

protected:
  F f;
  A a;
};

template < class F, class A >
struct Binder_3 {
  typedef typename F::result_type result_type;
  Binder_3(const F& f_, const A& a_) : f(f_), a(a_) {}

  template < class A1, class A2 >
  result_type
  operator()(const A1& a1, const A2& a2) const
  { return f(a1, a2, a); }

  template < class A1, class A2, class A3 >
  result_type
  operator()(const A1& a1, const A2& a2, const A3& a3) const
  { return f(a1, a2, a, a3); }

  template < class A1, class A2, class A3, class A4 >
  result_type
  operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
  { return f(a1, a2, a, a3, a4); }

protected:
  F f;
  A a;
};

template < class F, class A >
struct Binder_4 {
  typedef typename F::result_type result_type;
  Binder_4(const F& f_, const A& a_) : f(f_), a(a_) {}

  template < class A1, class A2, class A3 >
  result_type
  operator()(const A1& a1, const A2& a2, const A3& a3) const
  { return f(a1, a2, a3, a); }

  template < class A1, class A2, class A3, class A4 >
  result_type
  operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
  { return f(a1, a2, a3, a, a4); }

protected:
  F f;
  A a;
};

template < class F, class A >
struct Binder_5 {
  typedef typename F::result_type result_type;
  Binder_5(const F& f_, const A& a_) : f(f_), a(a_) {}

  template < class A1, class A2, class A3, class A4 >
  result_type
  operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
  { return f(a1, a2, a3, a4, a); }

protected:
  F f;
  A a;
};

namespace CGALi {
  template < int i >
  struct Binder_helper;

  template <>
  struct Binder_helper< 1 > {
    template < class T, class A >
    struct Help { typedef Binder_1< T, A > Type; };
  };

  template <>
  struct Binder_helper< 2 > {
    template < class T, class A >
    struct Help { typedef Binder_2< T, A > Type; };
  };

  template <>
  struct Binder_helper< 3 > {
    template < class T, class A >
    struct Help { typedef Binder_3< T, A > Type; };
  };

  template <>
  struct Binder_helper< 4 > {
    template < class T, class A >
    struct Help { typedef Binder_4< T, A > Type; };
  };

  template <>
  struct Binder_helper< 5 > {
    template < class T, class A >
    struct Help { typedef Binder_5< T, A > Type; };
  };
}

template < class T, class A, int i >
struct Binder_class {
#ifdef __sgi
  typedef typename CGALi::Binder_helper< i >::template
    Help< T, A >::Type Type;
#else
  typedef typename CGALi::Binder_helper< i >::Help< T, A >::Type Type;
#endif
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
  template < int i > struct N {};
}

template < class F1, class F2 >
struct Compose_1 {
  typedef typename F1::result_type result_type;
  enum { arity = F2::arity };

  Compose_1(const F1& f1_, const F2& f2_) : f1(f1_), f2(f2_) {}

  template < class A1 >
  result_type operator()() const
  { return f1(f2()); }

  template < class A1 >
  result_type operator()(const A1& a1) const
  { return f1(f2(a1)); }

  template < class A1, class A2 >
  result_type operator()(const A1& a1, const A2& a2) const
  { return f1(f2(a1, a2)); }

  template < class A1, class A2, class A3 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3) const
  { return f1(f2(a1, a2, a3)); }

  template < class A1, class A2, class A3, class A4 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4) const
  { return f1(f2(a1, a2, a3, a4)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4, const A5& a5) const
  { return f1(f2(a1, a2, a3, a4, a5)); }

protected:
  F1 f1;
  F2 f2;
};

template < class F1, class F2, class F3 >
struct Compose_2 {
  typedef typename F1::result_type result_type;
  enum { arity = F2::arity + F3::arity };

  Compose_2(const F1& f1_, const F2& f2_, const F3& f3_)
  : f1(f1_), f2(f2_), f3(f3_)
  {}

  result_type operator()() const
  { return call(CGALi::N<F2::arity>(),
                CGALi::N<F3::arity>(), CGALi::N<arity>()); }

  template < class A1 >
  result_type operator()(const A1& a1) const
  { return call(a1, CGALi::N<F2::arity>(),
                CGALi::N<F3::arity>(), CGALi::N<arity>()); }

  template < class A1, class A2 >
  result_type operator()(const A1& a1, const A2& a2) const
  { return call(a1, a2, CGALi::N<F2::arity>(),
                CGALi::N<F3::arity>(), CGALi::N<arity>()); }

  template < class A1, class A2, class A3 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3) const
  { return call(a1, a2, a3, CGALi::N<F2::arity>(),
                CGALi::N<F3::arity>(), CGALi::N<arity>()); }

  template < class A1, class A2, class A3, class A4 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4 ) const
  { return call(a1, a2, a3, a4, CGALi::N<F2::arity>(),
                CGALi::N<F3::arity>(), CGALi::N<arity>()); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4, const A5& a5) const
  { return call(a1, a2, a3, a4, a5, CGALi::N<F2::arity>(),
                CGALi::N<F3::arity>(), CGALi::N<arity>()); }


  result_type call(CGALi::N<1>, CGALi::N<0>, CGALi::N<1>) const
  { return f1(f2(), f3()); }

  template < class A1 >
  result_type call(const A1& a1,
                   CGALi::N<1>, CGALi::N<0>, CGALi::N<1>) const
  { return f1(f2(a1), f3()); }

  template < class A1 >
  result_type call(const A1& a1,
                   CGALi::N<0>, CGALi::N<1>, CGALi::N<1>) const
  { return f1(f2(), f3(a1)); }

  template < class A1, class A2 >
  result_type call(const A1& a1, const A2& a2,
                   CGALi::N<2>, CGALi::N<0>, CGALi::N<2>) const
  { return f1(f2(a1, a2), f3()); }

  template < class A1, class A2 >
  result_type call(const A1& a1, const A2& a2,
                   CGALi::N<1>, CGALi::N<1>, CGALi::N<2>) const
  { return f1(f2(a1), f3(a2)); }

  template < class A1, class A2 >
  result_type call(const A1& a1, const A2& a2,
                   CGALi::N<0>, CGALi::N<2>, CGALi::N<2>) const
  { return f1(f2(), f3(a1, a2)); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   CGALi::N<3>, CGALi::N<0>, CGALi::N<3>) const
  { return f1(f2(a1, a2, a3), f3()); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   CGALi::N<1>, CGALi::N<2>, CGALi::N<3>) const
  { return f1(f2(a1), f3(a2, a3)); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   CGALi::N<2>, CGALi::N<1>, CGALi::N<3>) const
  { return f1(f2(a1, a2), f3(a3)); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   CGALi::N<0>, CGALi::N<3>, CGALi::N<3>) const
  { return f1(f2(), f3(a1, a2, a3)); }

protected:
  F1 f1;
  F2 f2;
  F3 f3;
};

template < class F1, class F2 >
inline Compose_1< F1, F2 >
compose(const F1& f1, const F2& f2)
{ return Compose_1< F1, F2 >(f1, f2); }

template < class F1, class F2, class F3 >
inline Compose_2< F1, F2, F3 >
compose(const F1& f1, const F2& f2, const F3& f3)
{ return Compose_2< F1, F2, F3 >(f1, f2, f3); }

CGAL_END_NAMESPACE

#endif // CGAL_FUNCTIONAL_MSVC_H //
// EOF //
