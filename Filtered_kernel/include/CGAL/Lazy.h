// Copyright (c) 2005,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Sylvain Pion

#ifndef CGAL_LAZY_H
#define CGAL_LAZY_H

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/Handle.h>
#include <CGAL/Object.h>
#include <CGAL/Kernel/Type_mapper.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/min_max_n.h>
#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Default.h>
#include <CGAL/tss.h>
#include <CGAL/type_traits/is_iterator.h>
#include <CGAL/transforming_iterator.h>

#include <optional>
#include <variant>
#ifdef CGAL_LAZY_KERNEL_DEBUG
#  include <boost/optional/optional_io.hpp>
#endif

#include <boost/mpl/has_xxx.hpp>

#include <iostream>
#include <iterator>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <atomic>
#include <thread>
#include <mutex>

namespace CGAL {

template <typename D>
class Bbox_d;

template <class E,
          class A,
          class E2A,
          class K>
class Lazy_kernel_base;

template <typename AT, typename ET, typename E2A> class Lazy;

template <typename ET_>
class Lazy_exact_nt;

template <typename AT, typename ET, typename E2A>
inline
decltype(auto)
approx(const Lazy<AT,ET,E2A>& l)
{
  return l.approx();
}

template <typename AT, typename ET, typename E2A>
inline
const ET&
exact(const Lazy<AT,ET,E2A>& l)
{
  return l.exact();
}


template <typename AT, typename ET, typename E2A>
inline
int
depth(const Lazy<AT,ET,E2A>& l)
{
  return l.depth();
}


#define CGAL_LAZY_FORWARD(T) \
  inline const T & approx(const T& d) { return d; } \
  inline const T & exact (const T& d) { return d; } \
  inline int       depth (const T&  ) { return 0; }

CGAL_LAZY_FORWARD(Bbox_2)
CGAL_LAZY_FORWARD(Bbox_3)
#undef CGAL_LAZY_FORWARD

template<class A> Bbox_d<A> const& approx(Bbox_d<A> const& d) { return d; }
template<class A> Bbox_d<A> const& exact (Bbox_d<A> const& d) { return d; }
template<class A> int              depth (Bbox_d<A> const&  ) { return 0; }

template<class T>inline std::enable_if_t<std::is_arithmetic<T>::value||std::is_enum<T>::value, T> approx(T d){return d;}
template<class T>inline std::enable_if_t<std::is_arithmetic<T>::value||std::is_enum<T>::value, T> exact (T d){return d;}
template<class T>inline std::enable_if_t<std::is_arithmetic<T>::value||std::is_enum<T>::value, int> depth(T){return -1;}

template<class T>inline std::enable_if_t<std::is_arithmetic<T>::value, Quotient<T>> approx(Quotient<T> d){return d;}
template<class T>inline std::enable_if_t<std::is_arithmetic<T>::value, Quotient<T>> exact (Quotient<T> d){return d;}
template<class T>inline std::enable_if_t<std::is_arithmetic<T>::value, int> depth(Quotient<T>){return -1;}

// For tag classes: Return_base_tag, Homogeneous_tag, Null_vector, Origin
template<class T>inline std::enable_if_t<std::is_empty<T>::value, T> exact(T){return {};}
template<class T>inline std::enable_if_t<std::is_empty<T>::value, T> approx(T){return {};}
template<class T>inline std::enable_if_t<std::is_empty<T>::value, int> depth(T){return -1;}

namespace internal{

  template <typename AT, typename ET, typename E2A>
  struct Evaluate<Lazy<AT, ET, E2A>>
  {
      void operator()(const Lazy<AT, ET, E2A>& l)
      {
          exact(l);
      }
};

  template <typename ET>
  struct Evaluate<Lazy_exact_nt<ET>>
  {
      void operator()(const Lazy_exact_nt<ET>& l)
      {
          exact(l);
      }
  };
} // internal namespace

// For an iterator, exact/approx applies to the objects it points to
template <class T, class=std::enable_if_t<is_iterator_type<T,std::input_iterator_tag>::value>>
auto exact(T const& t) {return make_transforming_iterator(t,[](auto const&u)->decltype(auto){return CGAL::exact(u);});}
template <class T, class=std::enable_if_t<is_iterator_type<T,std::input_iterator_tag>::value>>
auto approx(T const& t) {return make_transforming_iterator(t,[](auto const&u)->decltype(auto){return CGAL::approx(u);});}
template <class T, class=std::enable_if_t<is_iterator_type<T,std::input_iterator_tag>::value>>
int depth(T const&) {return 1;} // FIXME: depth(*t) would be better when t is valid, but not for end iterators, and the true answer would iterate on the range, but we can't do that with only one iterator... We need to replace iterators with ranges to solve that.

#ifdef CGAL_LAZY_KERNEL_DEBUG
template <class T>
void
print_at(std::ostream& os, const T& at)
{
  os << at;
}

template <class T>
void
print_at(std::ostream& os, const std::vector<T>& at)
{
  os << "std::vector";
}

template <>
void
print_at(std::ostream& os, const Object& o)
{
  os << "Object";
}

template <class T1, class T2>
void
print_at(std::ostream& os, const std::pair<T1,T2> & at)
{
  os << "[ " << at.first << " | " << at.second << " ]" << std::endl ;
}


template <typename AT, typename ET, typename E2A>
inline
void
print_dag(const Lazy<AT,ET,E2A>& l, std::ostream& os, int level = 0)
{
  l.print_dag(os, level);
}

inline
void
print_dag(double d, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++)
    os << "    ";
  os << d << std::endl;
}

inline
void
msg(std::ostream& os, int level, const char* s)
{
    for(int i = 0; i < level; i++)
      os << "    ";
    os << s << std::endl;
}

inline
void
print_dag(const Null_vector&, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++)
    os << "    ";
  os << "Null_vector" << std::endl;
}

inline
void
print_dag(const Origin&, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++)
    os << "    ";
  os << "Origin" << std::endl;
}

inline
void
print_dag(const Return_base_tag&, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++)
    os << "    ";
  os << "Return_base_tag" << std::endl;
}
#endif


struct Depth_base {
#ifdef CGAL_PROFILE
  int depth_;

  Depth_base()
    : depth_(0)
  {}

  int depth() const { return depth_; }
  void set_depth(int i)
  {
    depth_ = i;
    CGAL_HISTOGRAM_PROFILER(std::string("[Lazy_kernel DAG depths]"), i);
                            //(unsigned) ::log2(double(i)));
  }
#else
  int depth() const { return 0; }
  void set_depth(int) {}
#endif
};

template<class T, bool=std::is_base_of<Handle, T>::value> struct Lazy_reset_member_1 {
  void operator()(T& t)const{ t = T(); }
};
template<class T> struct Lazy_reset_member_1<T, true> {
  void operator()(T& t)const{ t.reset(); }
};
template<class T>void lazy_reset_member(T&t) {
  Lazy_reset_member_1<T>()(t);
}
template<class...T, std::size_t...i>void lazy_reset_member_tuple(std::tuple<T...>&t, std::index_sequence<i...>) {
  auto ignore = [](auto&&...){};
  ignore ( (lazy_reset_member(std::get<i>(t)), 0) ... );
}
template<class...T>void lazy_reset_member(std::tuple<T...>&t) {
  lazy_reset_member_tuple(t, std::make_index_sequence<sizeof...(T)>());
}

// 0: safe default, AT is behind a pointer that can be atomically changed, and it doesn't disappear during update_exact
// 1: use plain AT without protection
// 2: split an interval as 2 atomic_double
// FIXME: CGAL_USE_SSE2 is clearly not the right condition
#ifdef CGAL_HAS_THREADS
template<class AT>struct Lazy_rep_selector { static constexpr int value = 0; };
# if defined CGAL_USE_SSE2 && !defined __SANITIZE_THREAD__ && !__has_feature(thread_sanitizer)
template<bool b>struct Lazy_rep_selector<Interval_nt<b>> { static constexpr int value = 1; };
template<bool b, int N>struct Lazy_rep_selector<std::array<Interval_nt<b>,N>> { static constexpr int value = 1; };
 // Need some declarations, including Simple_cartesian.h would also be possible.
 template<class>struct Simple_cartesian;
 template<class>class Point_2;
 template<class>class Point_3;
template<bool b>struct Lazy_rep_selector<CGAL::Point_2<CGAL::Simple_cartesian<CGAL::Interval_nt<b>>>> { static constexpr int value = 1; };
template<bool b>struct Lazy_rep_selector<CGAL::Point_3<CGAL::Simple_cartesian<CGAL::Interval_nt<b>>>> { static constexpr int value = 1; };
# else
template<bool b>struct Lazy_rep_selector<Interval_nt<b>> { static constexpr int value = 2; };
# endif
#else
template<class AT>struct Lazy_rep_selector { static constexpr int value = 1; };
#endif

template<class AT>
struct AT_wrap {
  AT at_;
  AT_wrap():at_(){}
  AT_wrap(AT const& a):at_(a){}
  AT_wrap(AT&& a):at_(std::move(a)){}
  AT const& at()const{return at_;}
};

// TODO: avoid initializing AT for nothing
template<class AT, class ET>
struct AT_ET_wrap : AT_wrap<AT> {
  ET et_;
  AT_ET_wrap():et_(){}
  AT_ET_wrap(ET const& e):et_(e){}
  AT_ET_wrap(ET&& e):et_(std::move(e)){}
  template<class A, class E>AT_ET_wrap(A&&a, E&&e):AT_wrap<AT>(std::forward<A>(a)),et_(std::forward<E>(e)){}
  ET const& et()const{return et_;}
};

// Abstract base class for lazy numbers and lazy objects
template <typename AT_, typename ET, typename E2A, int=Lazy_rep_selector<AT_>::value /* 0 */>
class Lazy_rep : public Rep, public Depth_base
{
  Lazy_rep (const Lazy_rep&) = delete; // cannot be copied.
  Lazy_rep& operator= (const Lazy_rep&) = delete; // cannot be copied.

public:

  typedef AT_ AT;
  typedef AT_ET_wrap<AT,ET> Indirect;

  AT_wrap<AT> at_orig{};
  mutable std::atomic<AT_wrap<AT>*> ptr_ { &at_orig };
  mutable std::once_flag once;

  Lazy_rep () {}

  template<class A>
  Lazy_rep (A&& a)
      : at_orig(std::forward<A>(a)){}

  template<class A>
  Lazy_rep (int count, A&& a)
    : Rep(count), at_orig(std::forward<A>(a)){}

  template<class A, class E>
  Lazy_rep (A&& a, E&& e)
      : ptr_(new AT_ET_wrap<AT,ET>(std::forward<A>(a), std::forward<E>(e))) {}

  AT const& approx() const
  {
    return ptr_.load(std::memory_order_consume)->at();
  }

  const ET & exact_unsafe() const
  {
    CGAL_assertion(!is_lazy());
    return static_cast<AT_ET_wrap<AT,ET>*>(ptr_.load(std::memory_order_relaxed))->et();
  }

  const ET & exact() const
  {
    // The test is unnecessary, only use it if benchmark says so, or in order to avoid calling Lazy_exact_Ex_Cst::update_exact() (which used to contain an assertion)
    //if (is_lazy())
    std::call_once(once, [this](){this->update_exact();});
    return exact_unsafe(); // call_once already synchronized memory
  }

  template<class A>
  void set_at(AT_ET_wrap<AT,ET>* p, A&& a) const {
    p->at_ = std::forward<A>(a);
  }
  void set_at(AT_ET_wrap<AT,ET>* p) const {
    p->at_ = E2A()(p->et());
  }
  void keep_at(AT_ET_wrap<AT,ET>* p) const {
    p->at_ = at_orig.at(); // do not move!
  }

  void set_ptr(AT_ET_wrap<AT,ET>* p) const {
    ptr_.store(p, std::memory_order_release);
  }

  // I think we should have different code for cases where there is some cleanup to do (say, a sum of 2 Lazy_exact_nt) and for cases where there isn't (a Lazy_exact_nt constructed from a double), but it may require making exact() virtual. Objects can be hidden in a tuple in Lazy_rep_n, so checking if there is something to clean requires some code. It isn't clear if we also need to restrict that to cases where update_exact doesn't touch AT. The special version would be basically: if(et==0){pet=new ET(...);if(!et.exchange(0,pet))delete pet; update at?}

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void print_at_et(std::ostream& os, int level) const
  {
    for(int i = 0; i < level; i++){
      os << "    ";
    }
    os << "Approximation: ";
    print_at(os, approx());
    os << std::endl;
    if(! is_lazy()){
      for(int i = 0; i < level; i++){
        os << "    ";
      }
      os << "Exact: ";
      print_at(os, exact_unsafe());
      os << std::endl;
#ifdef CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID
      for(int i = 0; i < level; i++){
        os << "    ";
      }
      os << "  (type: " << typeid(exact_unsafe()).name() << ")" << std::endl;
#endif // CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID
    }
  }

  virtual void print_dag(std::ostream& os, int level) const {}
#endif

  bool is_lazy() const { return ptr_.load(std::memory_order_relaxed) == &at_orig; }
  virtual void update_exact() const = 0;
  virtual ~Lazy_rep() {
#if !defined __SANITIZE_THREAD__ && !__has_feature(thread_sanitizer)
    auto* p = ptr_.load(std::memory_order_relaxed);
    if (p != &at_orig) {
      std::atomic_thread_fence(std::memory_order_acquire);
      delete static_cast<Indirect*>(p);
    }
#else
    auto* p = ptr_.load(std::memory_order_consume);
    if (p != &at_orig) delete static_cast<Indirect*>(p);
#endif
  }
};

/* How (un)safe is this? The goal is to minimize the overhead compared to a single-thread version by making the fast path almost identical.
 * For scalars on x86_64, the interval is aligned, so load/store instructions will not slice any double (although Intel does not explicitly guarantee it). On recent hardware, they should even be atomic, although without an official guarantee, and we don't need 128-bit atomicity anyway. The main danger is the unpredictable optimizations a compiler could apply (volatile would disable most of them, but it doesn't seem great), including replacing load/store with memcpy, where I fear some implementation/hardware combinations might slice sometimes. Making Interval_nt a pair of atomic_double would avoid this problem, but would likely incur a penalty since compilers don't optimize atomics much, and we shouldn't need to store/load all the time (TODO: benchmark).
 * For aggregate-like types (Simple_cartesian::Point_3), it should be ok for the same reason.
 * This is definitely NOT safe for a std::vector like a Point_d with Dynamic_dimension_tag, so it should only be enabled on a case by case basis, if at all. Storing a Point_3 piecewise with 6 atomic_double would be doable, but painful, and I didn't benchmark to check the performance. */
template <typename AT_, typename ET, typename E2A>
class Lazy_rep<AT_, ET, E2A, 1> : public Rep, public Depth_base
{
  Lazy_rep (const Lazy_rep&) = delete; // cannot be copied.
  Lazy_rep& operator= (const Lazy_rep&) = delete; // cannot be copied.

public:

  typedef AT_ AT;
  typedef ET Indirect;

  mutable AT at;
  mutable std::atomic<ET*> ptr_ { nullptr };
#ifdef CGAL_HAS_THREADS
  mutable std::once_flag once;
#endif

  Lazy_rep () {}

  template<class A>
  Lazy_rep (A&& a)
      : at(std::forward<A>(a)) {}

  template<class A>
  Lazy_rep (int count, A&& a)
    : Rep(count), at(std::forward<A>(a)){}

  template<class A, class E>
  Lazy_rep (A&& a, E&& e)
      : at(std::forward<A>(a)), ptr_(new ET(std::forward<E>(e))) {}

  AT const& approx() const
  {
    return at;
  }

  template<class A>
  void set_at(ET*, A&& a) const {
    at = std::forward<A>(a);
  }

  void set_at(ET* p) const {
    set_at(p, E2A()(*p));
  }
  void keep_at(ET*) const { }

  const ET & exact_unsafe() const
  {
    return *ptr_.load(std::memory_order_relaxed);
  }

  const ET & exact() const
  {
#ifdef CGAL_HAS_THREADS
    // The test is unnecessary, only use it if benchmark says so, or in order to avoid calling Lazy_exact_Ex_Cst::update_exact() (which used to contain an assertion)
    //if (is_lazy())
    std::call_once(once, [this](){this->update_exact();});
#else
    if (is_lazy())
      this->update_exact();
#endif
    return exact_unsafe(); // call_once already synchronized memory
  }

  void set_ptr(ET* p) const {
    ptr_.store(p, std::memory_order_release);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void print_at_et(std::ostream& os, int level) const
  {
    for(int i = 0; i < level; i++){
      os << "    ";
    }
    os << "Approximation: ";
    print_at(os, approx());
    os << std::endl;
    if(! is_lazy()){
      for(int i = 0; i < level; i++){
        os << "    ";
      }
      os << "Exact: ";
      print_at(os, exact_unsafe());
      os << std::endl;
#ifdef CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID
      for(int i = 0; i < level; i++){
        os << "    ";
      }
      os << "  (type: " << typeid(exact_unsafe()).name() << ")" << std::endl;
#endif // CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID
    }
  }

  virtual void print_dag(std::ostream& os, int level) const {}
#endif

  bool is_lazy() const { return ptr_.load(std::memory_order_relaxed) == nullptr; }
  virtual void update_exact() const = 0;
  virtual ~Lazy_rep() {
#if !defined __SANITIZE_THREAD__ && !__has_feature(thread_sanitizer)
    auto* p = ptr_.load(std::memory_order_relaxed);
    if (p != nullptr) {
      std::atomic_thread_fence(std::memory_order_acquire);
      delete p;
    }
#else
    auto* p = ptr_.load(std::memory_order_consume);
    if (p != nullptr) delete p;
#endif
  }
};

// do we need to (forward) declare Interval_nt?
template <bool b, typename ET, typename E2A>
class Lazy_rep<Interval_nt<b>, ET, E2A, 2> : public Rep, public Depth_base
{
  Lazy_rep (const Lazy_rep&) = delete; // cannot be copied.
  Lazy_rep& operator= (const Lazy_rep&) = delete; // cannot be copied.

public:

  typedef Interval_nt<b> AT;
  typedef ET Indirect;

  mutable std::atomic<double> x, y; // -inf, +sup
  mutable std::atomic<ET*> ptr_ { nullptr };
  mutable std::once_flag once;

  Lazy_rep (AT a = AT())
      : x(-a.inf()), y(a.sup()) {}

  template<class E>
  Lazy_rep (AT a, E&& e)
      : x(-a.inf()), y(a.sup()), ptr_(new ET(std::forward<E>(e))) {}

  AT approx() const
  {
    // Do not check that the interval is valid. Indeed, using IO/WKT/traits_point.h for instance,
    // one can default-construct a point, then set X, and then Y, which amounts to
    // Point_2(Point_2(X, Point_2().y()).x(), Y).
    // With Epeck, we have a default constructed array of Interval_nt in Point_2(),
    // then .y() returns a Lazy_exact_nt containing an invalid interval,
    // and when we read that interval we end up here.
    return AT(-x.load(std::memory_order_relaxed), y.load(std::memory_order_relaxed), typename AT::no_check_t());
  }

  void set_at(ET*, AT a) const {
    x.store(-a.inf(), std::memory_order_relaxed);
    y.store(a.sup(), std::memory_order_relaxed);
  }

  void set_at(ET* p) const {
    set_at(p, E2A()(*p));
  }
  void keep_at(ET*) const { }

  const ET & exact_unsafe() const
  {
    return *ptr_.load(std::memory_order_relaxed);
  }

  const ET & exact() const
  {
    // The test is unnecessary, only use it if benchmark says so, or in order to avoid calling Lazy_exact_Ex_Cst::update_exact() (which used to contain an assertion)
    //if (is_lazy())
    std::call_once(once, [this](){this->update_exact();});
    return exact_unsafe(); // call_once already synchronized memory
  }

  void set_ptr(ET* p) const {
    ptr_.store(p, std::memory_order_release);
  }

  // I think we should have different code for cases where there is some cleanup to do (say, a sum of 2 Lazy_exact_nt) and for cases where there isn't (a Lazy_exact_nt constructed from a double). Objects can be hidden in a tuple in Lazy_rep_n, so checking if there is something to clean requires some code. It isn't clear if we also need to restrict that to cases where update_exact doesn't touch AT. The special version would be basically: if(et==0){pet=new ET(...);if(!et.exchange(0,pet))delete pet; update at?}

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void print_at_et(std::ostream& os, int level) const
  {
    for(int i = 0; i < level; i++){
      os << "    ";
    }
    os << "Approximation: ";
    print_at(os, approx());
    os << std::endl;
    if(! is_lazy()){
      for(int i = 0; i < level; i++){
        os << "    ";
      }
      os << "Exact: ";
      print_at(os, exact_unsafe());
      os << std::endl;
#ifdef CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID
      for(int i = 0; i < level; i++){
        os << "    ";
      }
      os << "  (type: " << typeid(exact_unsafe()).name() << ")" << std::endl;
#endif // CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID
    }
  }

  virtual void print_dag(std::ostream& os, int level) const {}
#endif

  bool is_lazy() const { return ptr_.load(std::memory_order_relaxed) == nullptr; }
  virtual void update_exact() const = 0;
  virtual ~Lazy_rep() {
#if !defined __SANITIZE_THREAD__ && !__has_feature(thread_sanitizer)
    auto* p = ptr_.load(std::memory_order_relaxed);
    if (p != nullptr) {
      std::atomic_thread_fence(std::memory_order_acquire);
      delete p;
    }
#else
    auto* p = ptr_.load(std::memory_order_consume);
    if (p != nullptr) delete p;
#endif
  }
};


template<typename AT, typename ET, typename AC, typename EC, typename E2A, bool noprune, typename...L>
class Lazy_rep_n final :
  public Lazy_rep< AT, ET, E2A >, private EC
{
  typedef Lazy_rep< AT, ET, E2A > Base;
  // Lazy_rep_0 does not inherit from EC or take a parameter AC. It has different constructors.
  static_assert(sizeof...(L)>0, "Use Lazy_rep_0 instead");
  template <class Ei, class Ai, class E2Ai, class Ki> friend class Lazy_kernel_base;
  mutable std::tuple<L...> l; // L...l; is not yet allowed.
  const EC& ec() const { return *this; }
  template<std::size_t...I>
  void update_exact_helper(std::index_sequence<I...>) const {
    auto* p = new typename Base::Indirect(ec()( CGAL::exact( std::get<I>(l) ) ... ) );
    this->set_at(p);
    this->set_ptr(p);
    if(!noprune || is_currently_single_threaded())
      lazy_reset_member(l);
  }
  public:
  void update_exact() const {
    update_exact_helper(std::make_index_sequence<sizeof...(L)>{});
  }
  template<class...LL>
  Lazy_rep_n(const AC& ac, const EC& ec, LL&&...ll) :
    Lazy_rep<AT, ET, E2A>(ac(CGAL::approx(ll)...)), EC(ec), l(std::forward<LL>(ll)...)
  {
    this->set_depth((std::max)({ -1, (int)CGAL::depth(ll)...}) + 1);
  }
#ifdef CGAL_LAZY_KERNEL_DEBUG
  private:
  template<std::size_t...I>
  void print_dag_helper(std::ostream& os, int level, std::index_sequence<I...>) const {
    this->print_at_et(os, level);
    if(this->is_lazy()){
# ifdef CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID
      CGAL::msg(os, level, typeid(AC).name());
# endif
      CGAL::msg(os, level, "DAG with " "3" " child nodes:");
      using expander = int[];
      expander{0,(CGAL::print_dag(std::get<I>(l), os, level+1),0)...};
    }
  }
  public:
  void print_dag(std::ostream& os, int level) const {
    print_dag_helper(os, level, std::make_index_sequence<sizeof...(L)>{});
  }
#endif
};


template<typename AT, typename ET, typename AC, typename EC, typename E2A, typename...L>
class Lazy_rep_optional_n :
  public Lazy_rep< AT, ET, E2A >, private EC
{
  // Lazy_rep_0 does not inherit from EC or take a parameter AC. It has different constructors.
  static_assert(sizeof...(L)>0, "Use Lazy_rep_0 instead");
  template <class Ei, class Ai, class E2Ai, class Ki> friend class Lazy_kernel_base;
  mutable std::tuple<L...> l; // L...l; is not yet allowed.

  const EC& ec() const { return *this; }

  template<std::size_t...I>
  void update_exact_helper(std::index_sequence<I...>) const {
    typedef Lazy_rep< AT, ET, E2A > Base;
    auto* p = new typename Base::Indirect( * ec()( CGAL::exact( std::get<I>(l) ) ... ) );
    this->set_at(p);
    this->set_ptr(p);
    lazy_reset_member(l);
  }
  public:

  Lazy_rep_optional_n()
  {}

  void update_exact() const {
    update_exact_helper(std::make_index_sequence<sizeof...(L)>{});
  }

  template<class...LL>
  Lazy_rep_optional_n(const AT& a, const EC& ec, LL&&...ll) :
    Lazy_rep<AT, ET, E2A>(a), EC(ec), l(std::forward<LL>(ll)...)
  {
    this->set_depth((std::max)({ -1, (int)CGAL::depth(ll)...}) + 1);
  }

  template<class...LL>
  Lazy_rep_optional_n(int count, const AT& a, const EC& ec, LL&&...ll) :
    Lazy_rep<AT, ET, E2A>(count, a), EC(ec), l(std::forward<LL>(ll)...)
  {
    this->set_depth((std::max)({ -1, (int)CGAL::depth(ll)...}) + 1);
  }
#ifdef CGAL_LAZY_KERNEL_DEBUG
  private:
  template<std::size_t...I>
  void print_dag_helper(std::ostream& os, int level, std::index_sequence<I...>) const {
    this->print_at_et(os, level);
    if(this->is_lazy()){
# ifdef CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID
      CGAL::msg(os, level, typeid(AC).name());
# endif
      CGAL::msg(os, level, "DAG with " "3" " child nodes:");
      using expander = int[];
      expander{0,(CGAL::print_dag(std::get<I>(l), os, level+1),0)...};
    }
  }

  public:
  void print_dag(std::ostream& os, int level) const {
    print_dag_helper(os, level, std::make_index_sequence<sizeof...(L)>{});
  }
#endif
};

//____________________________________________________________
// The rep for the leaf node

template <typename AT, typename ET, typename E2A>
class Lazy_rep_0 final : public Lazy_rep<AT, ET, E2A>
{

  typedef Lazy_rep<AT, ET, E2A> Base;
public:

  void
  update_exact() const
  {
#ifdef CGAL_HAS_THREADS
    // Unless we add is_lazy before call_once in Lazy_rep. This test is
    // necessary because this class can be used either for default
    // construction, or to store a non-lazy exact value, and only the first one
    // should have a non-empty update_exact.
    // An alternative would be to add in the constructors taking an ET: std::call_once(this->once, [](){});
    if(!this->is_lazy()) return;
#endif
    auto* p = new typename Base::Indirect();
    this->set_ptr(p);
  }

  Lazy_rep_0()
    : Lazy_rep<AT,ET, E2A>() {}

  // TODO: the case where the exact value is provided at construction should
  // actually use a different class from the lazy default construction.
  template<class A, class E>
  Lazy_rep_0(A&& a, E&& e)
    : Lazy_rep<AT,ET,E2A>(std::forward<A>(a), std::forward<E>(e))
  {
    this->set_depth(0);
  }

#if 0
  // unused. Find a less ambiguous placeholder if necessary
  Lazy_rep_0(const AT& a, void*)
    : Lazy_rep<AT,ET,E2A>(a)
  {
    this->set_depth(0);
  }
#endif

  // E2A()(e) and std::forward<E>(e) could be evaluated in any order, but
  // that's ok, "forward" itself does not modify e, it may only mark it as
  // modifiable by the outer call, which is obviously sequenced after the inner
  // call E2A()(e).
  template<class E>
  Lazy_rep_0(E&& e)
    : Lazy_rep<AT,ET,E2A>(E2A()(e), std::forward<E>(e))
  {
    this->set_depth(0);
  }

  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
  }
};

#undef CGAL_LAZY_PRINT_TYPEID

template < typename K1, typename K2 >
struct Approx_converter
{
  typedef K1         Source_kernel;
  typedef K2         Target_kernel;
  //typedef Converter  Number_type_converter;

  template < typename T >
  decltype(auto)
  operator()(const T&t) const
  { return approx(t); }

  const Null_vector&
  operator()(const Null_vector& n) const
  { return n; }

  const Bbox_2&
  operator()(const Bbox_2& b) const
  { return b; }

  const Bbox_3&
  operator()(const Bbox_3& b) const
  { return b; }
};

template < typename K1, typename K2 >
struct Exact_converter
{
  typedef K1         Source_kernel;
  typedef K2         Target_kernel;
  //typedef Converter  Number_type_converter;

  template < typename T >
  decltype(auto)
  operator()(const T&t) const
  { return exact(t); }

  const Null_vector&
  operator()(const Null_vector& n) const
  { return n; }

  const Bbox_2&
  operator()(const Bbox_2& b) const
  { return b; }

  const Bbox_3&
  operator()(const Bbox_3& b) const
  { return b; }
};

//____________________________________________________________
// The handle class
template <typename AT_, typename ET_, typename E2A>
class Lazy : public Handle
{
  template <class Exact_kernel_,
            class Approximate_kernel_,
            class E2A_>
  friend struct Lazy_kernel;

  template <class E_,
            class A_,
            class E2A_,
            class K_>
  friend class Lazy_kernel_base;

public :

  typedef Lazy<AT_, ET_, E2A>  Self;
  typedef Lazy_rep<AT_, ET_, E2A>   Self_rep;

  typedef AT_ AT; // undocumented
  typedef ET_ ET; // undocumented

  typedef AT  Approximate_type;
  typedef ET  Exact_type;

/*
  typedef Self Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }
*/

  Lazy()
    : Handle(zero()) {}

  // Before Lazy::zero() used Boost.Thread, the definition of Lazy() was:
  //   Lazy()
  //   #ifndef CGAL_HAS_THREAD
  //     : Handle(zero()) {}
  //   #else
  //   {
  //     PTR = new Lazy_rep_0<AT, ET, E2A>();
  //   }
  //   #endif

  Lazy(Self_rep *r)
  {
    PTR = r;
  }

  Lazy(const ET& e)
  {
    PTR = new Lazy_rep_0<AT,ET,E2A>(e);
  }

  Lazy(ET&& e)
  {
    PTR = new Lazy_rep_0<AT,ET,E2A>(std::move(e));
  }

  friend void swap(Lazy& a, Lazy& b) noexcept
  { swap(static_cast<Handle&>(a), static_cast<Handle&>(b)); }

  decltype(auto) approx() const
  { return ptr()->approx(); }

  const ET& exact() const
  { return ptr()->exact(); }

  unsigned depth() const
  {
    return ptr()->depth();
  }

  void print_dag(std::ostream& os, int level) const
  {
    ptr()->print_dag(os, level);
  }

  private:

  // We have a static variable for optimizing the default constructor,
  // which is in particular heavily used for pruning DAGs.
  static const Self & zero()
  {
    // Note that the new only happens inside an if() inside the macro
    // So it would be a mistake to put the new before the macro
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(Self,z,(new Lazy_rep_0<AT, ET, E2A>()));
    return z;
  }

  Self_rep * ptr() const { return (Self_rep*) PTR; }
};

template <typename LK, typename AC, typename EC>
struct Lazy_construction_optional_for_polyhedral_envelope
{
  static const bool Protection = true;
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename LK::E2A E2A;
  typedef std::optional<typename LK::Point_3> result_type;

  CGAL_NO_UNIQUE_ADDRESS AC ac;
  CGAL_NO_UNIQUE_ADDRESS EC ec;


  // for Intersect_point_3 with 3 Plane_3
  template <typename L1>
  result_type operator()(const L1& l1, const L1& l2, const L1& l3) const
  {
    {
      Protect_FPU_rounding<Protection> P;
      try {
        std::optional<typename AK::Point_3> oap = ac(CGAL::approx(l1),CGAL::approx(l2),CGAL::approx(l3));
        if(oap == std::nullopt){
          return std::nullopt;
        }
        // Now we have to construct a rep for a lazy point with the three lazy planes.
        typedef Lazy_rep_optional_n<typename AK::Point_3, typename EK::Point_3, AC, EC, E2A, L1, L1, L1> LazyPointRep;
        CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(LazyPointRep, rep);

        const typename AK::Point_3 ap = *oap;
        // rep = LazyPointRep(2,ap, ec, l1, l2, l3);
        rep.~LazyPointRep(); new (&rep) LazyPointRep(2, ap, ec, l1, l2, l3);
        typename LK::Point_3 lp(&rep);
        return std::make_optional(lp);

      } catch (Uncertain_conversion_exception&) {}
    }
    Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
    std::optional<typename EK::Point_3> oep = ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3));
    if(oep == std::nullopt){
      return std::nullopt;
    }
    typedef Lazy_rep_0<typename AK::Point_3, typename EK::Point_3, E2A> LazyPointRep;
    const typename EK::Point_3 ep = *oep;
    LazyPointRep *rep = new LazyPointRep(ep);
    typename LK::Point_3 lp(rep);
    return std::make_optional(lp);
  }

  // for Intersect_point_3 with  Plane_3  Line_3
  template <typename L1, typename L2>
  result_type operator()(const L1& l1, const L2& l2) const
  {

    {
      Protect_FPU_rounding<Protection> P;
      try {
        std::optional<typename AK::Point_3> oap = ac(CGAL::approx(l1),CGAL::approx(l2));
        if(oap == std::nullopt){
          return std::nullopt;
        }
        // Now we have to construct a rep for a lazy point with the line and the plane.
        typedef Lazy_rep_optional_n<typename AK::Point_3, typename EK::Point_3, AC, EC, E2A, L1, L2> LazyPointRep;

        CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(LazyPointRep, rep);
        const typename AK::Point_3 ap = *oap;
        // rep = LazyPointRep(2, ap, ec, l1, l2);
        rep.~LazyPointRep(); new (&rep) LazyPointRep(2, ap, ec, l1, l2);
        typename LK::Point_3 lp(&rep);
        return std::make_optional(lp);

      } catch (Uncertain_conversion_exception&) {}
    }
    Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
    std::optional<typename EK::Point_3> oep = ec(CGAL::exact(l1), CGAL::exact(l2));
    if(oep == std::nullopt){
      return std::nullopt;
    }
    typedef Lazy_rep_0<typename AK::Point_3, typename EK::Point_3, E2A> LazyPointRep;
    const typename EK::Point_3 ep = *oep;
    LazyPointRep *rep = new LazyPointRep(ep);
    typename LK::Point_3 lp(rep);
    return std::make_optional(lp);
  }
};

// used in Newkernel_d
template <typename LK, typename AC, typename EC>
struct Lazy_construction_nt {
  Lazy_construction_nt(){}
  Lazy_construction_nt(LK const&){}

  static const bool Protection = true;

  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename LK::E2A E2A;

  CGAL_NO_UNIQUE_ADDRESS AC ac;
  CGAL_NO_UNIQUE_ADDRESS EC ec;

  template<class...L>
  auto operator()(L const&...l) const ->
  Lazy_exact_nt<CGAL::cpp20::remove_cvref_t<decltype(ec(CGAL::exact(l)...))>>
  {
    typedef CGAL::cpp20::remove_cvref_t<decltype(ec(CGAL::exact(l)...))> ET;
    typedef CGAL::cpp20::remove_cvref_t<decltype(ac(CGAL::approx(l)...))> AT;
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> P;
      try {
        return new Lazy_rep_n<AT, ET, AC, EC, To_interval<ET>, false, L... >(ac, ec, l...);
      } catch (Uncertain_conversion_exception&) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
    return new Lazy_rep_0<AT,ET,To_interval<ET> >(ec( CGAL::exact(l)... ));
  }
};


template <typename LK>
Object
make_lazy(const Object& eto)
{
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename LK::E2A E2A;

  if (eto.is_empty())
    return Object();

#define CGAL_Kernel_obj(X) \
  if (const typename EK::X* ptr = object_cast<typename EK::X>(&eto)) \
    return make_object(typename LK::X(new Lazy_rep_0<typename AK::X, typename EK::X, E2A>(*ptr)));

#include <CGAL/Kernel/interface_macros.h>

//now handle vector
#define CGAL_Kernel_obj(X) \
      {  \
        const std::vector<typename EK::X>* v_ptr;\
        if ( (v_ptr = object_cast<std::vector<typename EK::X> >(&eto)) ) { \
          std::vector<typename LK::X> V;\
          V.resize(v_ptr->size());                           \
          for (unsigned int i = 0; i < v_ptr->size(); ++i)                \
            V[i] = typename LK::X( new Lazy_rep_0<typename AK::X,typename EK::X,E2A>((*v_ptr)[i])); \
          return make_object(V);                                      \
        }\
      }

CGAL_Kernel_obj(Point_2)
CGAL_Kernel_obj(Point_3)
#undef CGAL_Kernel_obj


  std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's (#2)" << std::endl;
  std::cerr << "dynamic type of the Object : " << eto.type().name() << std::endl;

  return Object();
}

// This functor selects the i-th element in a vector of T2's
template <typename T2>
struct Ith_for_intersection {
  typedef T2 result_type;
  int i;

  Ith_for_intersection(int i_)
    : i(i_)
  {}

  const T2&
  operator()(const Object& o) const
  {
    const std::vector<T2>* ptr = object_cast<std::vector<T2> >(&o);
    return (*ptr)[i];
  }
};

template <typename LK, typename AC, typename EC>
struct Lazy_cartesian_const_iterator_2
{
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename LK::Cartesian_const_iterator_2 result_type;

  CGAL_NO_UNIQUE_ADDRESS AC ac;
  CGAL_NO_UNIQUE_ADDRESS EC ec;

public:

  template < typename L1>
  result_type
  operator()(const L1& l1) const
  {
    return result_type(&l1);
  }

  template < typename L1>
  result_type
  operator()(const L1& l1, int) const
  {
    return result_type(&l1,2);
  }

};


template <typename LK, typename AC, typename EC>
struct Lazy_cartesian_const_iterator_3
{
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename LK::Cartesian_const_iterator_3 result_type;

  CGAL_NO_UNIQUE_ADDRESS AC ac;
  CGAL_NO_UNIQUE_ADDRESS EC ec;

public:

  template < typename L1>
  result_type
  operator()(const L1& l1) const
  {
    return result_type(&l1);
  }

  template < typename L1>
  result_type
  operator()(const L1& l1, int) const
  {
    return result_type(&l1,3);
  }

};

template <typename T>
struct Object_cast
{
  typedef T result_type;

  const T&
  operator()(const Object& o) const
  {
    return *object_cast<T>(&o);
  }
};

// The following functor returns an Object with a Lazy<Something> inside
// As the nested kernels return Objects of AK::Something and EK::Something
// we have to unwrap them from the Object, and wrap them in a Lazy<Something>
template <typename LK, typename AC, typename EC>
struct Lazy_construction_object
{
  static const bool Protection = true;

  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename LK::E2A E2A;
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Object result_type;

  typedef Lazy<Object, Object, E2A> Lazy_object;

  CGAL_NO_UNIQUE_ADDRESS AC ac;
  CGAL_NO_UNIQUE_ADDRESS EC ec;

public:
  template <class... L>
  decltype(auto)
  operator()(const L&... l) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> P;
      try {
        Lazy_object lo(new Lazy_rep_n<result_type, result_type, AC, EC, E2A, false, L...>(ac, ec, l...));

        if(lo.approx().is_empty())
          return Object();

#define CGAL_Kernel_obj(X) \
        if (object_cast<typename AK::X>(& (lo.approx()))) { \
          typedef Lazy_rep_n< typename AK::X, typename EK::X, Object_cast<typename AK::X>, Object_cast<typename EK::X>, E2A, false, Lazy_object> Lcr; \
          Lcr * lcr = new Lcr(Object_cast<typename AK::X>(), Object_cast<typename EK::X>(), lo); \
          return make_object(typename LK::X(lcr)); \
        }

#include <CGAL/Kernel/interface_macros.h>

      // We now check vector<X>
#define CGAL_Kernel_obj(X) \
        { \
          const std::vector<typename AK::X>* v_ptr; \
          if ( (v_ptr = object_cast<std::vector<typename AK::X> >(& (lo.approx()))) ) { \
            std::vector<typename LK::X> V; \
            V.resize(v_ptr->size()); \
            for (unsigned int i = 0; i < v_ptr->size(); i++) { \
              V[i] = typename LK::X(new Lazy_rep_n<typename AK::X, typename EK::X, Ith_for_intersection<typename AK::X>, \
                                                   Ith_for_intersection<typename EK::X>, E2A, false, Lazy_object> \
                                    (Ith_for_intersection<typename AK::X>(i), Ith_for_intersection<typename EK::X>(i), lo)); \
            } \
            return make_object(V); \
          }\
        }

CGAL_Kernel_obj(Point_2)
CGAL_Kernel_obj(Point_3)
#undef CGAL_Kernel_obj

        std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's (#1)" << std::endl;
        std::cerr << "dynamic type of the Object : " << lo.approx().type().name() << std::endl;

        return Object();
      } catch (Uncertain_conversion_exception&) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
    ET eto = ec(CGAL::exact(l)...);
    return make_lazy<LK>(eto);
  }
};

//____________________________________________________________
// The magic functor that has Lazy<Something> as result type.

namespace internal {
BOOST_MPL_HAS_XXX_TRAIT_DEF(result_type)

// lift boost::get into a functor with a result_type member name and
// extend it to operate on optionals

// TODO there is a mismatch between the result_type typedef and the
// actual return type of operator()
template<typename T>
struct Variant_cast {
  typedef T result_type;

  template<typename ... U>
  const T&
  operator()(const std::optional< std::variant< U ... > >& o) const {
    // can throw but should never because we always build it inside
    // a static visitor with the right type
    return std::get<T>(*o);
  }

  template<typename ... U>
  T&
  operator()(std::optional< std::variant< U ... > >& o) const {
    // can throw but should never because we always build it inside
    // a static visitor with the right type, if it throws bad_get
    return std::get<T>(*o);
  }
};


template<typename Result, typename AK, typename LK, typename EK, typename Origin>
struct Fill_lazy_variant_visitor_2 {
  Fill_lazy_variant_visitor_2(Result& r, Origin& o) : r(&r), o(&o) {}
  Result* r;
  Origin* o;

  template<typename T>
  void operator()(const T&) {
    // the equivalent type we are currently matching in the lazy kernel
    typedef T AKT;
    typedef typename Type_mapper<AKT, AK, EK>::type EKT;
    typedef typename Type_mapper<AKT, AK, LK>::type LKT;

    typedef Lazy_rep_n<AKT, EKT, Variant_cast<AKT>, Variant_cast<EKT>, typename LK::E2A, false, Origin> Lcr;
    Lcr * lcr = new Lcr(Variant_cast<AKT>(), Variant_cast<EKT>(), *o);

    *r = LKT(lcr);
  }

  template<typename T>
  void operator()(const std::vector<T>& t) {
    typedef T AKT;
    typedef typename Type_mapper<AKT, AK, EK>::type EKT;
    typedef typename Type_mapper<AKT, AK, LK>::type LKT;

    std::vector<LKT> V;
    V.resize(t.size());
    for (unsigned int i = 0; i < t.size(); i++) {
      V[i] = LKT(new Lazy_rep_n<AKT, EKT, Ith_for_intersection<AKT>,
                 Ith_for_intersection<EKT>, typename LK::E2A, false, Origin>
                 (Ith_for_intersection<AKT>(i), Ith_for_intersection<EKT>(i), *o));
    }

    *r = V;
  }
};

template<typename Result, typename AK, typename LK, typename EK>
struct Fill_lazy_variant_visitor_0 {
  Fill_lazy_variant_visitor_0(Result& r) : r(&r) {}
  Result* r;

  template<typename T>
  void operator()(const T& t) {
    // the equivalent type we are currently matching in the lazy kernel
    typedef T EKT;
    typedef typename Type_mapper<EKT, EK, AK>::type AKT;
    typedef typename Type_mapper<EKT, EK, LK>::type LKT;

    *r = LKT(new Lazy_rep_0<AKT, EKT, typename LK::E2A>(t));
  }

  template<typename T>
  void operator()(const std::vector<T>& t) {
    typedef T EKT;
    typedef typename Type_mapper<EKT, EK, AK>::type AKT;
    typedef typename Type_mapper<EKT, EK, LK>::type LKT;

    std::vector<LKT> V;
    V.resize(t.size());
    for (unsigned int i = 0; i < t.size(); i++) {
      V[i] = LKT(new Lazy_rep_0<AKT, EKT, typename LK::E2A>(t[i]));
    }

    *r = V;
  }
};

} // internal

template<class AK, class AC>
struct Disable_lazy_pruning
{
  static const bool value = false;
};
template<class AK>
struct Disable_lazy_pruning<AK, typename AK::Construct_weighted_point_2>
{
  static const bool value = true;
};
template<class AK>
struct Disable_lazy_pruning<AK, typename AK::Construct_weighted_point_3>
{
  static const bool value = true;
};

template <typename LK, typename AC, typename EC, typename E2A_ = Default>
struct Lazy_construction
{
  static const bool Protection = true;

  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename Default::Get<E2A_, typename LK::E2A>::type E2A;

  CGAL_NO_UNIQUE_ADDRESS AC ac;
  CGAL_NO_UNIQUE_ADDRESS EC ec;

  // to detect the return type of Intersection_[23]'s functors
  template <typename T>
  struct is_optional_variant : std::false_type { };

  template <typename ...Args>
  struct is_optional_variant<std::optional<std::variant<Args...> > > : std::true_type { };

  template <class... L>
  decltype(auto)
  operator()(const L&... l) const
  {
    typedef typename Type_mapper<std::invoke_result_t<AC, typename Type_mapper<L, LK, AK>::type...>, AK, AK>::type AT;

    // ----------------------- FT -----------------------
    if constexpr (std::is_same_v<AT, typename AK::FT>)
    {
      typedef typename Type_mapper<std::invoke_result_t<EC,typename Type_mapper<L, LK, EK>::type...>,EK,EK>::type ET;

      typedef Lazy_exact_nt<CGAL::cpp20::remove_cvref_t<decltype(ec(CGAL::exact(l)...))>> result_type;

      CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
      {
        Protect_FPU_rounding<Protection> P;
        try {
          return result_type(new Lazy_rep_n<AT, ET, AC, EC, To_interval<ET>, false, L... >(ac, ec, l...));
        } catch (Uncertain_conversion_exception&) {}
      }
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
      return result_type(new Lazy_rep_0<AT,ET,To_interval<ET> >(ec( CGAL::exact(l)... )));
    }
    // ----------------------- Bounding boxes -----------------------
    else if constexpr (std::disjunction_v<std::is_same<AT, CGAL::Bbox_2>,
                                          std::is_same<AT, CGAL::Bbox_3> >)
    {
      CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
      {
        // Protection is outside the try block as VC8 has the CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG
        Protect_FPU_rounding<Protection> P;
        try {
          return ac(CGAL::approx(l...));
        } catch (Uncertain_conversion_exception&) {}
      }
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
      return ec(CGAL::exact(l...));
    }
    // ----------------------- CGAL::Object -----------------------
    else if constexpr (std::is_same_v<AT, CGAL::Object>)
    {
      typedef Lazy<Object, Object, E2A> Lazy_object;

      CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
      {
        Protect_FPU_rounding<Protection> P;
        try {
          Lazy_object lo(new Lazy_rep_n<CGAL::Object, CGAL::Object, AC, EC, E2A, false, L...>(ac, ec, l...));

          if(lo.approx().is_empty())
            return Object();

#         define CGAL_Kernel_obj(X) \
          if (object_cast<typename AK::X>(& (lo.approx()))) { \
            typedef Lazy_rep_n< typename AK::X, typename EK::X, Object_cast<typename AK::X>, Object_cast<typename EK::X>, E2A, false, Lazy_object> Lcr; \
            Lcr * lcr = new Lcr(Object_cast<typename AK::X>(), Object_cast<typename EK::X>(), lo); \
            return make_object(typename LK::X(lcr)); \
          }

#         include <CGAL/Kernel/interface_macros.h>

          // We now check vector<X>
#         define CGAL_Kernel_obj(X) \
          { \
            const std::vector<typename AK::X>* v_ptr; \
            if ( (v_ptr = object_cast<std::vector<typename AK::X> >(& (lo.approx()))) ) { \
              std::vector<typename LK::X> V; \
              V.resize(v_ptr->size()); \
              for (unsigned int i = 0; i < v_ptr->size(); i++) { \
                V[i] = typename LK::X(new Lazy_rep_n<typename AK::X, typename EK::X, Ith_for_intersection<typename AK::X>, \
                                                    Ith_for_intersection<typename EK::X>, E2A, false, Lazy_object> \
                                      (Ith_for_intersection<typename AK::X>(i), Ith_for_intersection<typename EK::X>(i), lo)); \
              } \
              return make_object(V); \
            }\
          }

          CGAL_Kernel_obj(Point_2)
          CGAL_Kernel_obj(Point_3)
#         undef CGAL_Kernel_obj

          std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's (#1)" << std::endl;
          std::cerr << "dynamic type of the Object : " << lo.approx().type().name() << std::endl;

          return Object();
        } catch (Uncertain_conversion_exception&) {}
      }
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
      CGAL::Object eto = ec(CGAL::exact(l)...);
      return make_lazy<LK>(eto);
    }
    // std::optional<std::variant< > > (Intersection_23 result types)
    else if constexpr (is_optional_variant<AT>::value)
    {
      typedef typename Type_mapper<std::invoke_result_t<EC, typename Type_mapper<L, LK, EK>::type...>,EK,EK>::type ET;

      typedef typename Type_mapper<std::invoke_result_t<AC, typename Type_mapper<L, LK, AK>::type...>, AK, LK>::type result_type;

      CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
      {
        Protect_FPU_rounding<Protection> P;

        try {
          Lazy<AT, ET, E2A> lazy(new Lazy_rep_n<AT, ET, AC, EC, E2A, false, L...>(AC(), EC(), l...));

          // the approximate result requires the trait with types from the AK
          AT approx_v = lazy.approx();
          // the result we build
          result_type res;

          if(!approx_v) {
            // empty
            return res;
          }

          // the static visitor fills the result_type with the correct unwrapped type
          internal::Fill_lazy_variant_visitor_2< result_type, AK, LK, EK, Lazy<AT, ET, E2A> > visitor(res, lazy);
          std::visit(visitor, *approx_v);

          return res;
        } catch (Uncertain_conversion_exception&) {}
      }
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
      ET exact_v = EC()(CGAL::exact(l)...);
      result_type res;

      if(!exact_v) {
        return res;
      }

      internal::Fill_lazy_variant_visitor_0<result_type, AK, LK, EK> visitor(res);
      std::visit(visitor, *exact_v);
      return res;
    }
    // ----------------------- GENERIC -----------------------
    else
    {
      typedef typename Type_mapper<std::invoke_result_t<EC,typename Type_mapper<L, LK, EK>::type...>,EK,EK>::type ET;

      typedef Lazy<AT, ET, E2A> Handle;
      typedef typename Type_mapper<std::invoke_result_t<AC, typename Type_mapper<L, LK, AK>::type...>, AK, LK>::type result_type;

      static const bool noprune = Disable_lazy_pruning<AK, AC>::value;

      CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
      {
        Protect_FPU_rounding<Protection> P;
        try {
          return result_type(Handle(new Lazy_rep_n<AT, ET, AC, EC, E2A, noprune, L...>(ac, ec, l...)));
        } catch (Uncertain_conversion_exception&) {}
      }
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
      return result_type(Handle(new Lazy_rep_0<AT, ET, E2A> (ec(CGAL::exact(l)...))));
    }
  }

  // nullary
  decltype(auto)
  operator()() const
  {
    typedef std::invoke_result_t<AC> AT;
    typedef std::invoke_result_t<EC> ET;
    typedef Lazy<AT, ET, E2A> Handle;
    typedef typename Type_mapper<AT, AK, LK>::type result_type;

    return result_type( Handle() );
  }
};

} //namespace CGAL


#include <CGAL/enable_warnings.h>

#endif // CGAL_LAZY_H
