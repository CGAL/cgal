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
// New Functor Adaptors (MSVC version).
// ============================================================================

#ifndef CGAL_FUNCTIONAL_MSVC_H
#define CGAL_FUNCTIONAL_MSVC_H 1

CGAL_BEGIN_NAMESPACE

// helper classes for arity
namespace CGALi {
#ifndef CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION

  template < class T > struct Arity_minus_one;
  template < int i > struct Arity_minus_one< Arity_tag< i > > {
    typedef Arity_tag< i - 1 > Arity;
  };

  template < class T1, class T2 > struct Arity_plus;
  template < int i, int j >
  struct Arity_plus< Arity_tag< i >, Arity_tag< j > > {
    typedef Arity_tag< i + j > Arity;
  };

#else
#ifdef CGAL_CFG_ENUM_BUG
#error "Too many compiler bugs, sorry ..."
#endif // CGAL_CFG_ENUM_BUG

  // msvc6 needs this, whysoever ...
  template < class T1, class T2 >
  struct CGAL__Wrap_arity {
    enum { a1 = T1::arity };
    enum { a2 = T2::arity };
    enum { ar1 = a1 - 1 };
    enum { ar2 = a1 + a2 };
  };

  template < class T > struct Arity_minus_one {
    typedef Arity_tag< CGAL__Wrap_arity< T, T >::ar1 > Arity;
  };

  template < class T1, class T2 > struct Arity_plus {
    typedef Arity_tag< CGAL__Wrap_arity< T1, T2 >::ar2 > Arity;
  };

#endif // CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION
} // namespace CGALi
namespace CGALi {

  template < class F >
  struct Swapper_1 {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper_1(const F& f_) : f(f_) {}
  
    template < class A1, class A2 >
    result_type operator()
    (const A1& a1, const A2& a2) const
    { return f(a2, a1); }
  
    template < class A1, class A2, class A3 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3) const
    { return f(a2, a1, a3); }
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return f(a2, a1, a3, a4); }
  
    template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3,
     const A4& a4, const A5& a5) const
    { return f(a2, a1, a3, a4, a5); }
  
  protected:
    F f;
  };
  
  template < class F >
  struct Swapper_2 {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper_2(const F& f_) : f(f_) {}
  
    template < class A1, class A2 >
    result_type operator()
    (const A1& a1, const A2& a2) const
    { return f(a1, a2); }
  
    template < class A1, class A2, class A3 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3) const
    { return f(a1, a3, a2); }
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return f(a1, a3, a2, a4); }
  
    template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3,
     const A4& a4, const A5& a5) const
    { return f(a1, a3, a2, a4, a5); }
  
  protected:
    F f;
  };
  
  template < class F >
  struct Swapper_3 {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper_3(const F& f_) : f(f_) {}
  
    template < class A1, class A2 >
    result_type operator()
    (const A1& a1, const A2& a2) const
    { return f(a1, a2); }
  
    template < class A1, class A2, class A3 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3) const
    { return f(a1, a2, a3); }
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return f(a1, a2, a4, a3); }
  
    template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3,
     const A4& a4, const A5& a5) const
    { return f(a1, a2, a4, a3, a5); }
  
  protected:
    F f;
  };
  
  template < class F >
  struct Swapper_4 {
    typedef typename F::result_type result_type;
    typedef typename F::Arity       Arity;
  
    Swapper_4(const F& f_) : f(f_) {}
  
    template < class A1, class A2 >
    result_type operator()
    (const A1& a1, const A2& a2) const
    { return f(a1, a2); }
  
    template < class A1, class A2, class A3 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3) const
    { return f(a1, a2, a3); }
  
    template < class A1, class A2, class A3, class A4 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return f(a1, a2, a3, a4); }
  
    template < class A1, class A2, class A3, class A4, class A5 >
    result_type operator()
    (const A1& a1, const A2& a2, const A3& a3,
     const A4& a4, const A5& a5) const
    { return f(a1, a2, a3, a5, a4); }
  
  protected:
    F f;
  };
  

  template < class F >
  struct Swap_helper {
    template < int i > struct Help;
    template <> struct Help< 1 > { typedef Swapper_1< F > Type; };
    template <> struct Help< 2 > { typedef Swapper_2< F > Type; };
    template <> struct Help< 3 > { typedef Swapper_3< F > Type; };
    template <> struct Help< 4 > { typedef Swapper_4< F > Type; };
  };

} // namespace CGALi

template < class F, int i = 1 >
struct Swap {
  typedef CGALi::Swap_helper< F >::Help< i >::Type Type;
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

template < class F, class A >
struct Binder_1 {
  typedef typename F::result_type result_type;
  typedef
    typename CGALi::Arity_minus_one< typename F::Arity >::Arity Arity;

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
  typedef
    typename CGALi::Arity_minus_one< typename F::Arity >::Arity Arity;

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
  typedef
    typename CGALi::Arity_minus_one< typename F::Arity >::Arity Arity;

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
  typedef
    typename CGALi::Arity_minus_one< typename F::Arity >::Arity Arity;

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
  typedef
    typename CGALi::Arity_minus_one< typename F::Arity >::Arity Arity;

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
  template < class T, class A >
  struct Binder_helper {
    template < int i > struct Help;
    template <> struct Help< 1 > { typedef Binder_1< T, A > Type; };
    template <> struct Help< 2 > { typedef Binder_2< T, A > Type; };
    template <> struct Help< 3 > { typedef Binder_3< T, A > Type; };
    template <> struct Help< 4 > { typedef Binder_4< T, A > Type; };
    template <> struct Help< 5 > { typedef Binder_5< T, A > Type; };
  };
}

template < class T, class A, int i >
struct Bind {
  typedef typename CGALi::Binder_helper< T, A >::Help< i >::Type Type;
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



template < class F1, class R1, class F2, class R2 >
struct Compose_1 {
  typedef typename F1::result_type result_type;
  typedef R2 Arity;

  Compose_1(const F1& f1_, const F2& f2_) : f1(f1_), f2(f2_) {}

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4, const A5& a5) const
  { return f1(f2(a1, a2, a3, a4, a5)); }

  template < class A1, class A2, class A3, class A4 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4) const
  { return f1(f2(a1, a2, a3, a4)); }

  template < class A1, class A2, class A3 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3) const
  { return f1(f2(a1, a2, a3)); }

  template < class A1, class A2 >
  result_type operator()(const A1& a1, const A2& a2) const
  { return f1(f2(a1, a2)); }

  result_type operator()() const
  { return f1(f2()); }

  template < class A1 >
  result_type operator()(const A1& a1) const
  { return f1(f2(a1)); }

protected:
  F1 f1;
  F2 f2;
};

template < class F1, class R1, class F2, class R2, class F3, class R3 >
struct Compose_2 {
  typedef typename F1::result_type result_type;
  typedef typename CGALi::Arity_plus< R2, R3 >::Arity Arity;

  Compose_2(const F1& f1_, const F2& f2_, const F3& f3_)
  : f1(f1_), f2(f2_), f3(f3_)
  {}

  result_type operator()() const
  { return call(R2(), R3(), Arity()); }

  template < class A1 >
  result_type operator()(const A1& a1) const
  { return call(a1, R2(), R3(), Arity()); }

  template < class A1, class A2 >
  result_type operator()(const A1& a1, const A2& a2) const
  { return call(a1, a2, R2(), R3(), Arity()); }

  template < class A1, class A2, class A3 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3) const
  { return call(a1, a2, a3, R2(), R3(), Arity()); }

  template < class A1, class A2, class A3, class A4 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4 ) const
  { return call(a1, a2, a3, a4, R2(), R3(),Arity()); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4, const A5& a5) const
  {
    return call(a1, a2, a3, a4, a5, R2(), R3(), Arity()); }

  result_type call(Arity_tag< 0 >, Arity_tag< 0 >, Arity_tag< 0 >) const
  { return f1(f2(), f3()); }

  template < class A1 >
  result_type call(const A1& a1,
                   Arity_tag< 1 >, Arity_tag< 0 >, Arity_tag< 1 >) const
  { return f1(f2(a1), f3()); }

  template < class A1 >
  result_type call(const A1& a1,
                   Arity_tag< 0 >, Arity_tag< 1 >, Arity_tag< 1 >) const
  { return f1(f2(), f3(a1)); }

  template < class A1, class A2 >
  result_type call(const A1& a1, const A2& a2,
                   Arity_tag< 2 >, Arity_tag< 0 >, Arity_tag< 2 >) const
  { return f1(f2(a1, a2), f3()); }

  template < class A1, class A2 >
  result_type call(const A1& a1, const A2& a2,
                   Arity_tag< 1 >, Arity_tag< 1 >, Arity_tag< 2 >) const
  { return f1(f2(a1), f3(a2)); }

  template < class A1, class A2 >
  result_type call(const A1& a1, const A2& a2,
                   Arity_tag< 0 >, Arity_tag< 2 >, Arity_tag< 2 >) const
  { return f1(f2(), f3(a1, a2)); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 3 >, Arity_tag< 0 >, Arity_tag< 3 >) const
  { return f1(f2(a1, a2, a3), f3()); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 2 >, Arity_tag< 1 >, Arity_tag< 3 >) const
  { return f1(f2(a1, a2), f3(a3)); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 1 >, Arity_tag< 2 >, Arity_tag< 3 >) const
  { return f1(f2(a1), f3(a2, a3)); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 0 >, Arity_tag< 3 >, Arity_tag< 3 >) const
  { return f1(f2(), f3(a1, a2, a3)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 4 >, Arity_tag< 0 >, Arity_tag< 4 >) const
  { return f1(f2(a1, a2, a3, a4), f3()); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 3 >, Arity_tag< 1 >, Arity_tag< 4 >) const
  { return f1(f2(a1, a2, a3), f3(a4)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 2 >, Arity_tag< 2 >, Arity_tag< 4 >) const
  { return f1(f2(a1, a2), f3(a3, a4)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 1 >, Arity_tag< 3 >, Arity_tag< 4 >) const
  { return f1(f2(a1), f3(a2, a3, a4)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 0 >, Arity_tag< 4 >, Arity_tag< 4 >) const
  { return f1(f2(), f3(a1, a2, a3, a4)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 5 >, Arity_tag< 0 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2, a3, a4, a5), f3()); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 4 >, Arity_tag< 1 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2, a3, a4), f3(a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 3 >, Arity_tag< 2 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2, a3), f3(a4, a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 2 >, Arity_tag< 3 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2), f3(a3, a4, a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 1 >, Arity_tag< 4 >, Arity_tag< 5 >) const
  { return f1(f2(a1), f3(a2, a3, a4, a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 0 >, Arity_tag< 5 >, Arity_tag< 5 >) const
  { return f1(f2(), f3(a1, a2, a3, a4, a5)); }

protected:
  F1 f1;
  F2 f2;
  F3 f3;
};

template < class F1, class F2, class F3 >
struct Compose_shared_2 {
  typedef typename F1::result_type result_type;
  typedef typename F2::Arity       Arity;

  Compose_shared_2(const F1& f1_, const F2& f2_, const F3& f3_)
  : f1(f1_), f2(f2_), f3(f3_)
  {}

  result_type operator()() const
  { return f1(f2(), f3()); }

  template < class A1 >
  result_type operator()(const A1& a1) const
  { return f1(f2(a1), f3(a1)); }

  template < class A1, class A2 >
  result_type operator()(const A1& a1, const A2& a2) const
  { return f1(f2(a1, a2), f3(a1, a2)); }

  template < class A1, class A2, class A3 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3) const
  { return f1(f2(a1, a2, a3), f3(a1, a2, a3)); }

  template < class A1, class A2, class A3, class A4 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4 ) const
  { return f1(f2(a1, a2, a3, a4), f3(a1, a2, a3, a4)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4, const A5& a5) const
  { return f1(f2(a1, a2, a3, a4, a5), f3(a1, a2, a3, a4, a5)); }

protected:
  F1 f1;
  F2 f2;
  F3 f3;
};

template < class F1, class R1, class F2, class R2,
           class F3, class R3, class F4, class R4 >
struct Compose_3 {
  typedef typename F1::result_type result_type;
  typedef typename CGALi::Arity_plus< R2, R3 >::Arity   Rtmp;
  typedef typename CGALi::Arity_plus< Rtmp, R4 >::Arity Arity;

  Compose_3(const F1& f1_, const F2& f2_, const F3& f3_,
            const F4& f4_)
  : f1(f1_), f2(f2_), f3(f3_), f4(f4_)
  {}

  result_type operator()() const
  { return call(R2(), R3(), R4(), Arity()); }

  template < class A1 >
  result_type operator()(const A1& a1) const
  { return call(a1, R2(), R3(), R4(), Arity()); }

  template < class A1, class A2 >
  result_type operator()(const A1& a1, const A2& a2) const
  { return call(a1, a2, R2(), R3(), R4(), Arity()); }

  template < class A1, class A2, class A3 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3) const
  { return call(a1, a2, a3, R2(), R3(), R4(), Arity()); }

  template < class A1, class A2, class A3, class A4 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4) const
  { return call(a1, a2, a3, a4, R2(), R3(), R4(), Arity()); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4, const A5& a5) const
  { return call(a1, a2, a3, a4, a5, R2(), R3(), R4(), Arity()); }

  result_type call(Arity_tag< 0 >, Arity_tag< 0 >,
                   Arity_tag< 0 >, Arity_tag< 0 >) const
  { return f1(f2(), f3(), f4()); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 5 >, Arity_tag< 0 >,
                   Arity_tag< 0 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2, a3, a4, a5), f3(), f4()); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 0 >, Arity_tag< 5 >,
                   Arity_tag< 0 >, Arity_tag< 5 >) const
  { return f1(f2(), f3(a1, a2, a3, a4, a5), f4()); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 0 >, Arity_tag< 0 >,
                   Arity_tag< 5 >, Arity_tag< 5 >) const
  { return f1(f2(), f3(), f4(a1, a2, a3, a4, a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 4 >, Arity_tag< 1 >,
                   Arity_tag< 0 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2, a3, a4), f3(a5), f4()); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 4 >, Arity_tag< 0 >,
                   Arity_tag< 1 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2, a3, a4), f3(), f4(a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 1 >, Arity_tag< 4 >,
                   Arity_tag< 0 >, Arity_tag< 5 >) const
  { return f1(f2(a1), f3(a2, a3, a4, a5), f4()); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 0 >, Arity_tag< 4 >,
                   Arity_tag< 1 >, Arity_tag< 5 >) const
  { return f1(f2(), f3(a1, a2, a3, a4), f4(a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 1 >, Arity_tag< 0 >,
                   Arity_tag< 4 >, Arity_tag< 5 >) const
  { return f1(f2(a1), f3(), f4(a2, a3, a4, a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 0 >, Arity_tag< 1 >,
                   Arity_tag< 4 >, Arity_tag< 5 >) const
  { return f1(f2(), f3(a1), f4(a2, a3, a4, a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 3 >, Arity_tag< 2 >,
                   Arity_tag< 0 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2, a3), f3(a4, a5), f4()); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 3 >, Arity_tag< 0 >,
                   Arity_tag< 2 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2, a3), f3(), f4(a4, a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 2 >, Arity_tag< 3 >,
                   Arity_tag< 0 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2), f3(a3, a4, a5), f4()); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 0 >, Arity_tag< 3 >,
                   Arity_tag< 2 >, Arity_tag< 5 >) const
  { return f1(f2(), f3(a1, a2, a3), f4(a4, a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 2 >, Arity_tag< 0 >,
                   Arity_tag< 3 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2), f3(), f4(a3, a4, a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 0 >, Arity_tag< 2 >,
                   Arity_tag< 3 >, Arity_tag< 5 >) const
  { return f1(f2(), f3(a1, a2), f4(a3, a4, a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 3 >, Arity_tag< 1 >,
                   Arity_tag< 1 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2, a3), f3(a4), f4(a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 1 >, Arity_tag< 3 >,
                   Arity_tag< 1 >, Arity_tag< 5 >) const
  { return f1(f2(a1), f3(a2, a3, a4), f4(a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 1 >, Arity_tag< 1 >,
                   Arity_tag< 3 >, Arity_tag< 5 >) const
  { return f1(f2(a1), f3(a2), f4(a3, a4, a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 2 >, Arity_tag< 2 >,
                   Arity_tag< 1 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2), f3(a3, a4), f4(a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 2 >, Arity_tag< 1 >,
                   Arity_tag< 2 >, Arity_tag< 5 >) const
  { return f1(f2(a1, a2), f3(a3), f4(a4, a5)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4, const A5& a5,
                   Arity_tag< 1 >, Arity_tag< 2 >,
                   Arity_tag< 2 >, Arity_tag< 5 >) const
  { return f1(f2(a1), f3(a2, a3), f4(a4, a5)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 4 >, Arity_tag< 0 >,
                   Arity_tag< 0 >, Arity_tag< 4 >) const
  { return f1(f2(a1, a2, a3, a4), f3(), f4()); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 0 >, Arity_tag< 4 >,
                   Arity_tag< 0 >, Arity_tag< 4 >) const
  { return f1(f2(), f3(a1, a2, a3, a4), f4()); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 0 >, Arity_tag< 0 >,
                   Arity_tag< 4 >, Arity_tag< 4 >) const
  { return f1(f2(), f3(), f4(a1, a2, a3, a4)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 3 >, Arity_tag< 1 >,
                   Arity_tag< 0 >, Arity_tag< 4 >) const
  { return f1(f2(a1, a2, a3), f3(a4), f4()); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 3 >, Arity_tag< 0 >,
                   Arity_tag< 1 >, Arity_tag< 4 >) const
  { return f1(f2(a1, a2, a3), f3(), f4(a4)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 1 >, Arity_tag< 3 >,
                   Arity_tag< 0 >, Arity_tag< 4 >) const
  { return f1(f2(a1), f3(a2, a3, a4), f4()); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 0 >, Arity_tag< 3 >,
                   Arity_tag< 1 >, Arity_tag< 4 >) const
  { return f1(f2(), f3(a1, a2, a3), f4(a4)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 1 >, Arity_tag< 0 >,
                   Arity_tag< 3 >, Arity_tag< 4 >) const
  { return f1(f2(a1), f3(), f4(a2, a3, a4)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 0 >, Arity_tag< 1 >,
                   Arity_tag< 3 >, Arity_tag< 4 >) const
  { return f1(f2(), f3(a1), f4(a2, a3, a4)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 2 >, Arity_tag< 2 >,
                   Arity_tag< 0 >, Arity_tag< 4 >) const
  { return f1(f2(a1, a2), f3(a3, a4), f4()); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 2 >, Arity_tag< 0 >,
                   Arity_tag< 2 >, Arity_tag< 4 >) const
  { return f1(f2(a1, a2), f3(), f4(a3, a4)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 0 >, Arity_tag< 2 >,
                   Arity_tag< 2 >, Arity_tag< 4 >) const
  { return f1(f2(), f3(a1, a2), f4(a3, a4)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 2 >, Arity_tag< 1 >,
                   Arity_tag< 1 >, Arity_tag< 4 >) const
  { return f1(f2(a1, a2), f3(a3), f4(a4)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 1 >, Arity_tag< 2 >,
                   Arity_tag< 1 >, Arity_tag< 4 >) const
  { return f1(f2(a1), f3(a2, a3), f4(a4)); }

  template < class A1, class A2, class A3, class A4 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   const A4& a4,
                   Arity_tag< 1 >, Arity_tag< 1 >,
                   Arity_tag< 2 >, Arity_tag< 4 >) const
  { return f1(f2(a1), f3(a2), f4(a3, a4)); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 3 >, Arity_tag< 0 >,
                   Arity_tag< 0 >, Arity_tag< 3 >) const
  { return f1(f2(a1, a2, a3), f3(), f4()); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 0 >, Arity_tag< 3 >,
                   Arity_tag< 0 >, Arity_tag< 3 >) const
  { return f1(f2(), f3(a1, a2, a3), f4()); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 0 >, Arity_tag< 0 >,
                   Arity_tag< 3 >, Arity_tag< 3 >) const
  { return f1(f2(), f3(), f4(a1, a2, a3)); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 2 >, Arity_tag< 1 >,
                   Arity_tag< 0 >, Arity_tag< 3 >) const
  { return f1(f2(a1, a2), f3(a3), f4()); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 2 >, Arity_tag< 0 >,
                   Arity_tag< 1 >, Arity_tag< 3 >) const
  { return f1(f2(a1, a2), f3(), f4(a3)); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 1 >, Arity_tag< 2 >,
                   Arity_tag< 0 >, Arity_tag< 3 >) const
  { return f1(f2(a1), f3(a2, a3), f4()); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 0 >, Arity_tag< 2 >,
                   Arity_tag< 1 >, Arity_tag< 3 >) const
  { return f1(f2(), f3(a1, a2), f4(a3)); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 1 >, Arity_tag< 0 >,
                   Arity_tag< 2 >, Arity_tag< 3 >) const
  { return f1(f2(a1), f3(), f4(a2, a3)); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 0 >, Arity_tag< 1 >,
                   Arity_tag< 2 >, Arity_tag< 3 >) const
  { return f1(f2(), f3(a1), f4(a2, a3)); }

  template < class A1, class A2, class A3 >
  result_type call(const A1& a1, const A2& a2, const A3& a3,
                   Arity_tag< 1 >, Arity_tag< 1 >,
                   Arity_tag< 1 >, Arity_tag< 3 >) const
  { return f1(f2(a1), f3(a2), f4(a3)); }

  template < class A1, class A2 >
  result_type call(const A1& a1, const A2& a2,
                   Arity_tag< 2 >, Arity_tag< 0 >,
                   Arity_tag< 0 >, Arity_tag< 2 >) const
  { return f1(f2(a1, a2), f3(), f4()); }

  template < class A1, class A2 >
  result_type call(const A1& a1, const A2& a2,
                   Arity_tag< 0 >, Arity_tag< 2 >,
                   Arity_tag< 0 >, Arity_tag< 2 >) const
  { return f1(f2(), f3(a1, a2), f4()); }

  template < class A1, class A2 >
  result_type call(const A1& a1, const A2& a2,
                   Arity_tag< 0 >, Arity_tag< 0 >,
                   Arity_tag< 2 >, Arity_tag< 2 >) const
  { return f1(f2(), f3(), f4(a1, a2)); }

  template < class A1, class A2 >
  result_type call(const A1& a1, const A2& a2,
                   Arity_tag< 1 >, Arity_tag< 1 >,
                   Arity_tag< 0 >, Arity_tag< 2 >) const
  { return f1(f2(a1), f3(a2), f4()); }

  template < class A1, class A2 >
  result_type call(const A1& a1, const A2& a2,
                   Arity_tag< 1 >, Arity_tag< 0 >,
                   Arity_tag< 1 >, Arity_tag< 2 >) const
  { return f1(f2(a1), f3(), f4(a2)); }

  template < class A1, class A2 >
  result_type call(const A1& a1, const A2& a2,
                   Arity_tag< 0 >, Arity_tag< 1 >,
                   Arity_tag< 1 >, Arity_tag< 2 >) const
  { return f1(f2(), f3(a1), f4(a2)); }

  template < class A1 >
  result_type call(const A1& a1,
                   Arity_tag< 1 >, Arity_tag< 0 >,
                   Arity_tag< 0 >, Arity_tag< 1 >) const
  { return f1(f2(a1), f3(), f4()); }

  template < class A1 >
  result_type call(const A1& a1,
                   Arity_tag< 0 >, Arity_tag< 1 >,
                   Arity_tag< 0 >, Arity_tag< 1 >) const
  { return f1(f2(), f3(a1), f4()); }

  template < class A1 >
  result_type call(const A1& a1,
                   Arity_tag< 0 >, Arity_tag< 0 >,
                   Arity_tag< 1 >, Arity_tag< 1 >) const
  { return f1(f2(), f3(), f4(a1)); }

protected:
  F1 f1;
  F2 f2;
  F3 f3;
  F4 f4;
};

template < class F1, class F2, class F3, class F4 >
struct Compose_shared_3 {
  typedef typename F1::result_type result_type;
  typedef typename F2::Arity       Arity;

  Compose_shared_3(const F1& f1_, const F2& f2_,
                   const F3& f3_, const F4& f4_)
  : f1(f1_), f2(f2_), f3(f3_), f4(f4_)
  {}

  result_type operator()() const
  { return f1(f2(), f3(), f4()); }

  template < class A1 >
  result_type operator()(const A1& a1) const
  { return f1(f2(a1), f3(a1), f4(a1)); }

  template < class A1, class A2 >
  result_type operator()(const A1& a1, const A2& a2) const
  { return f1(f2(a1, a2), f3(a1, a2), f4(a1, a2)); }

  template < class A1, class A2, class A3 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3) const
  { return f1(f2(a1, a2, a3), f3(a1, a2, a3), f4(a1, a2, a3)); }

  template < class A1, class A2, class A3, class A4 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4 ) const
  { return f1(f2(a1, a2, a3, a4),
              f3(a1, a2, a3, a4),
              f4(a1, a2, a3, a4)); }

  template < class A1, class A2, class A3, class A4, class A5 >
  result_type operator()(const A1& a1, const A2& a2, const A3& a3,
                         const A4& a4, const A5& a5) const
  { return f1(f2(a1, a2, a3, a4, a5),
              f3(a1, a2, a3, a4, a5),
              f4(a1, a2, a3, a4, a5)); }

protected:
  F1 f1;
  F2 f2;
  F3 f3;
  F4 f4;
};


// ------------------------------------------------------------------------
// Encapsulate the composition type ==> can be adapted
// ------------------------------------------------------------------------

namespace CGALi {
  // msvc6 needs this, whysoever ...
  template < class T >
  struct CGAL__Wrap_type {
    typedef typename T::Arity Arity;
  };

  template < class F1, class F2, class F3, class F4 >
  struct Compose_helper {
    template < class T > struct Help;

    template <> struct Help< Arity_tag< 1 > >
    { typedef Compose_1< F1, Arity_tag< 1 >,
                         F2, typename CGAL__Wrap_type< F2 >::Arity
                     > Type; };

    template <> struct Help< Arity_tag< 2 > >
    { typedef Compose_2< F1, Arity_tag< 2 >,
                         F2, typename CGAL__Wrap_type< F2 >::Arity,
                         F3, typename CGAL__Wrap_type< F3 >::Arity
                     > Type; };

    template <> struct Help< Arity_tag< 3 > >
    { typedef Compose_3< F1, Arity_tag< 3 >,
                         F2, typename CGAL__Wrap_type< F2 >::Arity,
                         F3, typename CGAL__Wrap_type< F3 >::Arity,
                         F4, typename CGAL__Wrap_type< F4 >::Arity
                     > Type; };
  };

  template < class F1, class F2, class F3, class F4 >
  struct Compose_shared_helper {
    template < class T > struct Help;

    template <> struct Help< Arity_tag< 2 > >
    { typedef Compose_shared_2< F1, F2, F3 > Type; };

    template <> struct Help< Arity_tag< 3 > >
    { typedef Compose_shared_3< F1, F2, F3, F4 > Type; };
  };
}

template < class F1, class F2, class F3 = F2, class F4 = F2 >
struct Compose {
  typedef typename CGALi::Compose_helper< F1, F2, F3, F4 >::Help<
    typename CGALi::CGAL__Wrap_type< F1 >::Arity  >::Type Type;
};

template < class F1, class F2, class F3 = F2, class F4 = F2 >
struct Compose_shared {
  typedef typename CGALi::Compose_shared_helper< F1, F2, F3, F4 >::Help<
    typename CGALi::CGAL__Wrap_type< F1 >::Arity  >::Type Type;
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

#endif // CGAL_FUNCTIONAL_MSVC_H //
// EOF //
