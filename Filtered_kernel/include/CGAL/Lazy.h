// Copyright (c) 2005,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Andreas Fabri, Sylvain Pion

#ifndef CGAL_LAZY_H
#define CGAL_LAZY_H

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
#include <vector>
#include <CGAL/Default.h>

#ifdef CGAL_HAS_THREADS
#  include <boost/thread/tss.hpp>
#endif

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <boost/mpl/has_xxx.hpp>

#include <boost/preprocessor/facilities/expand.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/enum.hpp>

namespace CGAL {

template <typename AT, typename ET, typename EFT, typename E2A> class Lazy;

template <typename ET_>
class Lazy_exact_nt;

template <typename AT, typename ET, typename EFT, typename E2A>
inline
const AT&
approx(const Lazy<AT,ET, EFT, E2A>& l)
{
  return l.approx();
}

// Where is this one (non-const) needed ?  Is it ?
template <typename AT, typename ET, typename EFT, typename E2A>
inline
AT&
approx(Lazy<AT,ET, EFT, E2A>& l)
{
  return l.approx();
}


template <typename AT, typename ET, typename EFT, typename E2A>
inline
const ET&
exact(const Lazy<AT,ET,EFT,E2A>& l)
{
  return l.exact();
}


template <typename AT, typename ET, typename EFT, typename E2A>
inline
unsigned
depth(const Lazy<AT,ET,EFT,E2A>& l)
{
  return l.depth();
}


#define CGAL_LAZY_FORWARD(T) \
  inline const T & approx(const T& d) { return d; } \
  inline const T & exact (const T& d) { return d; } \
  inline unsigned  depth (const T&  ) { return 0; }


CGAL_LAZY_FORWARD(long double)
CGAL_LAZY_FORWARD(double)
CGAL_LAZY_FORWARD(float)
CGAL_LAZY_FORWARD(int)
CGAL_LAZY_FORWARD(unsigned int)
CGAL_LAZY_FORWARD(long)
CGAL_LAZY_FORWARD(unsigned long)
#ifdef CGAL_USE_LONG_LONG
CGAL_LAZY_FORWARD(long long)
CGAL_LAZY_FORWARD(unsigned long long)
#endif
CGAL_LAZY_FORWARD(Return_base_tag)
CGAL_LAZY_FORWARD(Null_vector)
CGAL_LAZY_FORWARD(Origin)
CGAL_LAZY_FORWARD(Orientation)
CGAL_LAZY_FORWARD(Bbox_2)
CGAL_LAZY_FORWARD(Bbox_3)

#undef CGAL_LAZY_FORWARD

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


template <typename AT, typename ET, typename EFT, typename E2A>
inline
void
print_dag(const Lazy<AT,ET,EFT,E2A>& l, std::ostream& os, int level = 0)
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
  unsigned depth_;
  Depth_base() { set_depth(0); }
  unsigned depth() const { return depth_; }
  void set_depth(unsigned i)
  {
    depth_ = i;
    CGAL_HISTOGRAM_PROFILER(std::string("[Lazy_kernel DAG depths]"), i);
                            //(unsigned) ::log2(double(i)));
  }
#else
  unsigned depth() const { return 0; }
  void set_depth(unsigned) {}
#endif
};


// Abstract base class for lazy numbers and lazy objects
template <typename AT_, typename ET, typename E2A>
class Lazy_rep : public Rep, public Depth_base
{
  Lazy_rep (const Lazy_rep&); // cannot be copied.

public:

  typedef AT_ AT;

  mutable AT at;
  mutable ET *et;

  Lazy_rep ()
    : at(), et(NULL){}

  Lazy_rep (const AT& a)
      : at(a), et(NULL){}

  Lazy_rep (const AT& a, const ET& e)
      : at(a), et(new ET(e)) {}

  const AT& approx() const
  {
      return at;
  }

  AT& approx()
  {
      return at;
  }

  const ET & exact() const
  {
    if (et==NULL)
      update_exact();
    return *et;
  }

  ET & exact()
  {
    if (et==NULL)
      update_exact();
    return *et;
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void print_at_et(std::ostream& os, int level) const
  {
    for(int i = 0; i < level; i++){
      os << "    ";
    }
    os << "Approximation: ";
    print_at(os, at);
    os << std::endl;
    if(! is_lazy()){
      for(int i = 0; i < level; i++){
	os << "    ";
      }
      os << "Exact: ";
      print_at(os, *et);
      os << std::endl;
#ifdef CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID
      for(int i = 0; i < level; i++){
	os << "    ";
      }
      os << "  (type: " << typeid(*et).name() << ")" << std::endl;
#endif // CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID
    }
  }

  virtual void print_dag(std::ostream& os, int level) const {}
#endif

  bool is_lazy() const { return et == NULL; }
  virtual void update_exact() const = 0;
  virtual ~Lazy_rep() { delete et; }
};


//____________________________________________________________
// The rep for the leaf node

template <typename AT, typename ET, typename E2A>
class Lazy_rep_0 : public Lazy_rep<AT, ET, E2A>
{

  typedef Lazy_rep<AT, ET, E2A> Base;
public:

  void
  update_exact() const
  {
    this->et = new ET();
  }

  Lazy_rep_0()
    : Lazy_rep<AT,ET, E2A>() {}

  Lazy_rep_0(const AT& a, const ET& e)
    : Lazy_rep<AT,ET,E2A>(a, e) {}

  Lazy_rep_0(const AT& a, void*)
    : Lazy_rep<AT,ET,E2A>(a) {}

  Lazy_rep_0(const ET& e)
    : Lazy_rep<AT,ET,E2A>(E2A()(e), e) {}

  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
  }
};

// Macro helpers to build the kernel objects
#define CGAL_TYPEMAP_AC(z, n, t) typedef typename Type_mapper< t##n, LK, AK >::type A##n;
#define CGAL_TYPEMAP_EC(z, n, t) typedef typename Type_mapper< t##n, LK, EK >::type E##n;
#define CGAL_LEXACT(z,n,t) CGAL::exact( l##n )
#define CGAL_LARGS(z, n, t) L##n const& l##n

#define CGAL_TMAP(z, n, d) typename Type_mapper< L##n, d##K, LK >::type
#define CGAL_PRUNE_TREE(z, n, d) l##n = L##n ();
#define CGAL_LINIT(z, n, d) l##n(l##n)
#define CGAL_LN(z, n, d) d(l##n)
#define CGAL_MLIST(z, n, d) mutable L##n l##n;

//____________________________________________________________

template <typename AT, typename ET, typename AC, typename EC, typename E2A, typename L1>
class Lazy_rep_1
  : public Lazy_rep<AT, ET, E2A>
  , private EC
{
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET(ec()(CGAL::exact(l1_)));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
  }

  Lazy_rep_1(const AC& ac, const EC& ec, const L1& l1)
    : Lazy_rep<AT,ET, E2A>(ac(CGAL::approx(l1))), EC(ec), l1_(l1)
  {
    this->set_depth(CGAL::depth(l1_) + 1);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID
#  define CGAL_LAZY_PRINT_TYPEID CGAL::msg(os, level, typeid(AC).name());
#else  // not CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID
#  define CGAL_LAZY_PRINT_TYPEID
#endif // not CGAL_LAZY_KERNEL_DEBUG_SHOW_TYPEID

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      CGAL_LAZY_PRINT_TYPEID
      CGAL::msg(os, level, "DAG with one child node:");
      CGAL::print_dag(l1_, os, level+1);
    }
  }
#endif
};

#ifdef CGAL_LAZY_KERNEL_DEBUG
#  define CGAL_PRINT_DAG_LN(z, n, d) \
  CGAL::print_dag(l##n, os, level+1);
#  define CGAL_LAZY_REP_PRINT_DAG(n)                            \
  void print_dag(std::ostream& os, int level) const {           \
    this->print_at_et(os, level);                               \
    if(this->is_lazy()){                                        \
      CGAL_LAZY_PRINT_TYPEID                                    \
      CGAL::msg(os, level, "DAG with " #n " child nodes:");     \
      BOOST_PP_REPEAT(n, CGAL_PRINT_DAG_LN, _)                  \
    }                                                           \
  }
#else // not CGAL_LAZY_KERNEL_DEBUG
#  define CGAL_LAZY_REP_PRINT_DAG(n)
#endif // not CGAL_LAZY_KERNEL_DEBUG

#define CGAL_LAZY_REP(z, n, d)                                               \
  template< typename AT, typename ET, typename AC, typename EC, typename E2A, BOOST_PP_ENUM_PARAMS(n, typename L)> \
class Lazy_rep_##n :public Lazy_rep< AT, \
                                     ET, \
                                     E2A >,                             \
                    private EC                                          \
{                                                                       \
  BOOST_PP_REPEAT(n, CGAL_MLIST, _)                                          \
  const EC& ec() const { return *this; } \
public: \
  void update_exact() const { \
    this->et = new ET(ec()( BOOST_PP_ENUM(n, CGAL_LEXACT, _) ) ); \
    this->at = E2A()(*(this->et));                           \
    BOOST_PP_REPEAT(n, CGAL_PRUNE_TREE, _) \
  } \
  Lazy_rep_##n(const AC& ac, const EC&, BOOST_PP_ENUM(n, CGAL_LARGS, _)) \
    : Lazy_rep<AT, ET, E2A>(ac( BOOST_PP_ENUM(n, CGAL_LN, CGAL::approx) )), BOOST_PP_ENUM(n, CGAL_LINIT, _) \
  { this->set_depth(max_n( BOOST_PP_ENUM(n, CGAL_LN, CGAL::depth) ) + 1); }  \
                                                                        \
  CGAL_LAZY_REP_PRINT_DAG(n)                                          \
};

BOOST_PP_REPEAT_FROM_TO(2, 9, CGAL_LAZY_REP, _)

#undef CGAL_TMAP
#undef CGAL_PRUNE_TREE
#undef CGAL_LINIT
#undef CGAL_LAZY_REP
#undef CGAL_LN
#undef CGAL_MLIST
#undef CGAL_PRINT_DAG_LN
#undef CGAL_LAZY_REP_PRINT_DAG
#undef CGAL_LAZY_PRINT_TYPEID

template < typename K1, typename K2 >
struct Approx_converter
{
  typedef K1         Source_kernel;
  typedef K2         Target_kernel;
  //typedef Converter  Number_type_converter;

  template < typename T >
  const typename T::AT&
  operator()(const T&t) const
  { return t.approx(); }

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
  const typename T::ET&
  operator()(const T&t) const
  { return t.exact(); }

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



template <typename AC, typename EC, typename E2A, typename L1>
class Lazy_rep_with_vector_1
  : public Lazy_rep<std::vector<Object>, std::vector<Object>, E2A>
  , private EC
{
  typedef std::vector<Object> AT;
  typedef std::vector<Object> ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const 
  {
// TODO : This looks really unfinished...
    std::vector<Object> vec;
    this->et = new ET();
    //this->et->reserve(this->at.size());
    ec()(CGAL::exact(l1_), std::back_inserter(*(this->et)));
    if(this->et==NULL)
    E2A()(*(this->et));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
  }

  Lazy_rep_with_vector_1(const AC& ac, const EC& /*ec*/, const L1& l1)
    : l1_(l1)
  {
    ac(CGAL::approx(l1), std::back_inserter(this->at));
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    os << "A Lazy_rep_with_vector_1 of size " <<  this->at.size() << std::endl;
    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with one child node:");
      CGAL::print_dag(l1_, os, level+1);

    }
  }
#endif
};


template <typename AC, typename EC, typename E2A, typename L1, typename L2>
class Lazy_rep_with_vector_2
  : public Lazy_rep<std::vector<Object>, std::vector<Object>, E2A>
  , private EC
{
  typedef std::vector<Object> AT;
  typedef std::vector<Object> ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET();
    this->et->reserve(this->at.size());
    ec()(CGAL::exact(l1_), CGAL::exact(l2_), std::back_inserter(*(this->et)));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
  }

  Lazy_rep_with_vector_2(const AC& ac, const EC& /*ec*/, const L1& l1, const L2& l2)
    : l1_(l1), l2_(l2)
  {
    ac(CGAL::approx(l1), CGAL::approx(l2), std::back_inserter(this->at));
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    os << "A Lazy_rep_with_vector_2 of size " <<  this->at.size() << std::endl;
    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with two child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
    }
  }
#endif
};


template <typename AC, typename EC, typename E2A, typename L1, typename L2, typename R1>
class Lazy_rep_2_1
  : public Lazy_rep<typename R1::AT, typename R1::ET, E2A>
  , private EC
{
  typedef typename R1::AT AT;
  typedef typename R1::ET ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET();
    ec()(CGAL::exact(l1_), CGAL::exact(l2_), *(this->et));
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
  }

  Lazy_rep_2_1(const AC& ac, const EC& /*ec*/, const L1& l1, const L2& l2)
    : Lazy_rep<AT,ET,E2A>(), l1_(l1), l2_(l2)
  {
    ac(CGAL::approx(l1), CGAL::approx(l2), this->at);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    os << "A Lazy_rep_2_1" << std::endl;
    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with two child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
    }
  }
#endif
};


//____________________________________________________________________________________
// The following rep class stores two non-const reference parameters of type R1 and R2

template <typename AC, typename EC, typename E2A, typename L1, typename L2, typename R1, typename R2>
class Lazy_rep_2_2
  : public Lazy_rep<std::pair<typename R1::AT,typename R2::AT>, std::pair<typename R1::ET, typename R2::ET>, E2A>
  , private EC
{
  typedef std::pair<typename R1::AT, typename R2::AT> AT;
  typedef std::pair<typename R1::ET, typename R2::ET> ET;
  typedef Lazy_rep<AT, ET, E2A> Base;

  mutable L1 l1_;
  mutable L2 l2_;

  const EC& ec() const { return *this; }

public:

  void
  update_exact() const
  {
    this->et = new ET();
    ec()(CGAL::exact(l1_), CGAL::exact(l2_), this->et->first, this->et->second );
    this->at = E2A()(*(this->et));
    // Prune lazy tree
    l1_ = L1();
    l2_ = L2();
  }

  Lazy_rep_2_2(const AC& ac, const EC& /*ec*/, const L1& l1, const L2& l2)
    : Lazy_rep<AT,ET,E2A>(), l1_(l1), l2_(l2)
  {
    ac(CGAL::approx(l1), CGAL::approx(l2), this->at.first, this->at.second);
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    os << "A Lazy_rep_2_2"  << std::endl;
    if(this->is_lazy()){
      CGAL::msg(os, level, "DAG with two child nodes:");
      CGAL::print_dag(l1_, os, level+1);
      CGAL::print_dag(l2_, os, level+1);
    }
  }
#endif
};


//____________________________________________________________
// The handle class
template <typename AT_, typename ET_, typename EFT, typename E2A>
class Lazy : public Handle
{
public :

  typedef Lazy<AT_, ET_, EFT, E2A>  Self;
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

  const AT& approx() const
  { return ptr()->approx(); }

  const ET& exact() const
  { return ptr()->exact(); }

  AT& approx()
  { return ptr()->approx(); }

  ET& exact()
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
#ifdef CGAL_HAS_THREADS
    static boost::thread_specific_ptr<Self> z;
    if (z.get() == NULL) {
        z.reset(new Self(new Lazy_rep_0<AT, ET, E2A>()));
    }
    return * z.get();
#else
    static const Self z = new Lazy_rep_0<AT, ET, E2A>();
    return z;
#endif
  }

  Self_rep * ptr() const { return (Self_rep*) PTR; }
};

// The magic functor for Construct_bbox_[2,3], as there is no Lazy<Bbox>

template <typename LK, typename AC, typename EC>
struct Lazy_construction_bbox
{
  static const bool Protection = true;
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename AC::result_type result_type;

  AC ac;
  EC ec;

  template <typename L1>
  result_type operator()(const L1& l1) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    // Protection is outside the try block as VC8 has the CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG
    Protect_FPU_rounding<Protection> P;
    try {
      return ac(CGAL::approx(l1));
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      return ec(CGAL::exact(l1));
    }
  }
};


template <typename LK, typename AC, typename EC>
struct Lazy_construction_nt {

  static const bool Protection = true;

  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename LK::E2A E2A;

  AC ac;
  EC ec;

  template<typename>
  struct result { };

#define CGAL_RESULT_NT(z, n, d)                                              \
  template< typename F, BOOST_PP_ENUM_PARAMS(n, class T) >              \
  struct result<F( BOOST_PP_ENUM_PARAMS(n, T) )> {                      \
    BOOST_PP_REPEAT(n, CGAL_TYPEMAP_EC, T)                                   \
    typedef Lazy_exact_nt<                                              \
      typename boost::remove_cv< typename boost::remove_reference <     \
      typename cpp11::result_of<EC( BOOST_PP_ENUM_PARAMS(n, E) )>::type >::type >::type > type; \
  };

  BOOST_PP_REPEAT_FROM_TO(1, 6, CGAL_RESULT_NT, _)

#define CGAL_NT_OPERATOR(z, n, d)                                            \
  template<BOOST_PP_ENUM_PARAMS(n, class L)>                            \
  typename cpp11::result_of<Lazy_construction_nt(BOOST_PP_ENUM_PARAMS(n, L))>::type \
  operator()( BOOST_PP_ENUM(n, CGAL_LARGS, _) ) const {                      \
    BOOST_PP_REPEAT(n, CGAL_TYPEMAP_EC, L)                                     \
    BOOST_PP_REPEAT(n, CGAL_TYPEMAP_AC, L)                                     \
    typedef typename boost::remove_cv< typename boost::remove_reference < \
                                        typename cpp11::result_of< EC(BOOST_PP_ENUM_PARAMS(n, E)) >::type >::type >::type ET; \
    typedef typename boost::remove_cv< typename boost::remove_reference < \
                                        typename cpp11::result_of< AC(BOOST_PP_ENUM_PARAMS(n, A)) >::type >::type >::type AT; \
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp); \
    Protect_FPU_rounding<Protection> P;                                 \
    try {                                                               \
      return new Lazy_rep_##n<AT, ET, AC, EC, To_interval<ET>, BOOST_PP_ENUM_PARAMS(n, L) >(ac, ec, BOOST_PP_ENUM_PARAMS(n, l)); \
    } catch (Uncertain_conversion_exception) {                          \
      CGAL_BRANCH_PROFILER_BRANCH(tmp);                                 \
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);          \
      return new Lazy_rep_0<AT,ET,To_interval<ET> >(ec( BOOST_PP_ENUM(n, CGAL_LEXACT, _) )); \
    }                                                                   \
  }                                                                     \

  BOOST_PP_REPEAT_FROM_TO(1, 6, CGAL_NT_OPERATOR, _)

#undef INTERVAL_OPERATOR
#undef CGAL_RESULT_NT
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


// This functor selects the i'th element in a vector of Object's
// and casts it to what is in the Object

template <typename T2>
struct Ith {
  typedef T2 result_type;

  // We keep a Sign member object
  // for future utilisation, in case
  // we have pairs of 2 T2 objects e.g.
  // for a numeric_point vector returned
  // from a construction of a possible
  // lazy algebraic kernel

  int i;
  Sign sgn;

  Ith(int i_)
    : i(i_)
  {sgn=NEGATIVE;}

  Ith(int i_, bool b_)
    : i(i_)
  { sgn= (b_) ? POSITIVE : ZERO;}

  const T2&
  operator()(const std::vector<Object>& v) const
  {
    if(sgn==NEGATIVE)
    return *object_cast<T2>(&v[i]);

    typedef std::pair<T2,unsigned int >         Pair_type_1;
    typedef std::pair<T2,std::pair<bool,bool> > Pair_type_2;

    if(const Pair_type_1 *p1 = object_cast<Pair_type_1>(&v[i]))
    	return p1->first;
    else if(const Pair_type_2 *p2 = object_cast<Pair_type_2>(&v[i]))
        return p2->first;

    CGAL_error_msg( " Unexpected encapsulated type ");
  }
};

// This functor selects the i'th element in a vector of T2's
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

// This functor selects the i'th element in a vector of T2's
template <typename T2>
struct Ith_for_intersection_with_variant {
  typedef T2 result_type;
  int i;

  Ith_for_intersection_with_variant(int i_)
    : i(i_)
  {}

  template< BOOST_VARIANT_ENUM_PARAMS(typename U) >
  const T2&
  operator()(const boost::optional< boost::variant< BOOST_VARIANT_ENUM_PARAMS(U) > >& o) const
  {
    const std::vector<T2>* ptr = (boost::get<std::vector<T2> >(&(*o)));
    return (*ptr)[i];
  }

  template< BOOST_VARIANT_ENUM_PARAMS(typename U) >
  const T2&
  operator()(const boost::variant< BOOST_VARIANT_ENUM_PARAMS(U) >& o) const
  {
    const std::vector<T2>* ptr = (boost::get<std::vector<T2> >(&o));
    return (*ptr)[i];
  }
};

template <typename LK, typename AC, typename EC>
struct Lazy_cartesian_const_iterator_2
{
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename LK::Cartesian_const_iterator_2 result_type;

  AC ac;
  EC ec;

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

  AC ac;
  EC ec;

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


// This is the magic functor for functors that write their result in a  reference argument
// In a first version we assume that the references are of type Lazy<Something>,
// and that the result type is void

template <typename LK, typename AK, typename EK, typename AC, typename EC, typename EFT, typename E2A>
struct Lazy_functor_2_1
{
  static const bool Protection = true;
  typedef void result_type;

  AC ac;
  EC ec;

public:

  template <typename L1, typename L2, typename R1>
  void
  operator()(const L1& l1, const L2& l2, R1& r1) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      // we suppose that R1 is a Lazy<Something>
      r1 = R1(new Lazy_rep_2_1<AC, EC, E2A, L1, L2, R1>(ac, ec, l1, l2));
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      typename R1::ET et;
      ec(CGAL::exact(l1), CGAL::exact(l2), et);
      r1 = R1(new Lazy_rep_0<typename R1::AT,typename R1::ET,E2A>(et));
    }
  }
};


template <typename T>
struct First
{
   typedef typename T::first_type result_type;

   const typename T::first_type&
   operator()(const T& p) const
   {
     return p.first;
   }
 };

template <typename T>
struct Second
{
  typedef typename T::second_type result_type;

  const typename T::second_type&
  operator()(const T& p) const
  {
    return p.second;
  }
};

// This is the magic functor for functors that write their result in a reference argument
// In a first version we assume that the references are of type Lazy<Something>,
// and that the result type is void

//template <typename LK, typename AK, typename EK, typename AC, typename EC, typename EFT, typename E2A>
template <typename LK, typename AC, typename EC>
struct Lazy_functor_2_2
{
  static const bool Protection = true;

  typedef void result_type;
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename EK::FT EFT;
  typedef typename LK::E2A E2A;

  AC ac;
  EC ec;

public:

  template <typename L1, typename L2, typename R1, typename R2>
  void
  operator()(const L1& l1, const L2& l2, R1& r1, R2& r2) const
  {
    typedef Lazy<typename R1::AT, typename R1::ET, EFT, E2A> Handle_1;
    typedef Lazy<typename R2::AT, typename R2::ET, EFT, E2A> Handle_2;
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      typedef Lazy<std::pair<typename R1::AT, typename R2::AT>, std::pair<typename R1::ET, typename R2::ET>, EFT, E2A> Lazy_pair;
      Lazy_pair lv(new Lazy_rep_2_2<AC, EC, E2A, L1, L2, R1, R2>(ac, ec, l1, l2));
      // lv->approx() is a std::pair<R1::AT, R2::AT>;
      r1 = R1(Handle_1(new Lazy_rep_1<void, void, First<std::pair<typename R1::AT, typename R2::AT> >, First<std::pair<typename R1::ET, typename R2::ET> >, E2A, Lazy_pair>(First<std::pair<typename R1::AT, typename R2::AT> >(), First<std::pair<typename R1::ET, typename R2::ET> >(), lv)));
      r2 = R2(Handle_2(new Lazy_rep_1<void, void, Second<std::pair<typename R1::AT, typename R2::AT> >, Second<std::pair<typename R1::ET, typename R2::ET> >, E2A, Lazy_pair>(Second<std::pair<typename R1::AT, typename R2::AT> >(), Second<std::pair<typename R1::ET, typename R2::ET> >(), lv)));
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      typename R1::ET et1, et2;
      ec(CGAL::exact(l1), CGAL::exact(l2), et1, et2);
      r1 = R1(Handle_1(new Lazy_rep_0<typename R1::AT,typename R1::ET,E2A>(et1)));
      r2 = R2(Handle_2(new Lazy_rep_0<typename R2::AT,typename R2::ET,E2A>(et2)));
    }
  }
};


// This is the magic functor for functors that write their result as Objects into an output iterator

template <typename LK, typename AC, typename EC>
struct Lazy_intersect_with_iterators
{
  static const bool Protection = true;
  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename EK::FT EFT;
  typedef typename LK::E2A E2A;
  typedef void result_type;
  typedef Lazy<Object, Object, EFT, E2A> Lazy_object;
  typedef Lazy<std::vector<Object>, std::vector<Object>, EFT, E2A> Lazy_vector;

  AC ac;
  EC ec;

public:

  // In the example we intersect two Lazy<Segment>s
  // and write into a back_inserter(list<Object([Lazy<Point>,Lazy<Segment>]) >)
  template <typename L1, typename L2, typename OutputIterator>
  OutputIterator
  operator()(const L1& l1, const L2& l2, OutputIterator it) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      Lazy_vector lv(new Lazy_rep_with_vector_2<AC, EC, E2A, L1, L2>(ac, ec, l1, l2));
      // lv.approx() is a std::vector<Object([AK::Point_2,AK::Segment_2])>
      // that is, when we get here we have constructed all approximate results
      for (unsigned int i = 0; i < lv.approx().size(); i++) {
// FIXME : I'm not sure how this work...
#define CGAL_Kernel_obj(X) if (object_cast<typename AK::X>(& (lv.approx()[i]))) { \
	  *it++ = make_object(typename LK::X(new Lazy_rep_1<typename AK::X, typename EK::X, Ith<typename AK::X>, \
                                                                      Ith<typename EK::X>, E2A, Lazy_vector> \
                                                 (Ith<typename AK::X>(i), Ith<typename EK::X>(i), lv))); \
          continue; \
	}

#include <CGAL/Kernel/interface_macros.h>

        std::cerr << "we need  more casts" << std::endl;
      }

    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      // TODO: Instead of using a vector, write an iterator adapter
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      std::vector<Object> exact_objects;
      ec(CGAL::exact(l1), CGAL::exact(l2), std::back_inserter(exact_objects));
      for (std::vector<Object>::const_iterator oit = exact_objects.begin();
	   oit != exact_objects.end();
	   ++oit){
	*it++ = make_lazy<LK>(*oit);
      }
    }
    return it;
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
//
// TODO: write operators for other than two arguments. For the current kernel we only need two for Intersect_2

template <typename LK, typename AC, typename EC>
struct Lazy_construction_object
{
  static const bool Protection = true;

  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename EK::FT EFT;
  typedef typename LK::E2A E2A;
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Object result_type;

  typedef Lazy<Object, Object, EFT, E2A> Lazy_object;

  AC ac;
  EC ec;

public:

  template <typename L1>
  result_type
  operator()(const L1& l1) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      Lazy_object lo(new Lazy_rep_1<result_type, result_type, AC, EC, E2A, L1>(ac, ec, l1));

      if(lo.approx().is_empty())
        return Object();

#define CGAL_Kernel_obj(X) \
      if (object_cast<typename AK::X>(& (lo.approx()))) { \
	typedef Lazy_rep_1< typename AK::X, typename EK::X, Object_cast<typename AK::X>, Object_cast<typename EK::X>, E2A, Lazy_object> Lcr; \
	Lcr * lcr = new Lcr(Object_cast<typename AK::X>(), Object_cast<typename EK::X>(), lo); \
	return make_object(typename LK::X(lcr)); \
      }

#include <CGAL/Kernel/interface_macros.h>

      std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's (#1)" << std::endl;
      std::cerr << "dynamic type of the Object : " << lo.approx().type().name() << std::endl;

    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      ET eto = ec(CGAL::exact(l1));
      return make_lazy<LK>(eto);
    }
    return Object();
  }

  template <typename L1, typename L2>
  result_type
  operator()(const L1& l1, const L2& l2) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      Lazy_object lo(new Lazy_rep_2<result_type, result_type, AC, EC, E2A, L1, L2>(ac, ec, l1, l2));

      if(lo.approx().is_empty())
        return Object();

#define CGAL_Kernel_obj(X) \
      if (object_cast<typename AK::X>(& (lo.approx()))) { \
	typedef Lazy_rep_1<typename AK::X, typename EK::X, Object_cast<typename AK::X>, Object_cast<typename EK::X>, E2A, Lazy_object> Lcr; \
	Lcr * lcr = new Lcr(Object_cast<typename AK::X>(), Object_cast<typename EK::X>(), lo); \
	return make_object(typename LK::X(lcr)); \
      }

#include <CGAL/Kernel/interface_macros.h>

    // We now check vector<X>
  
#define CGAL_Kernel_obj(X) \
      {  \
        const std::vector<typename AK::X>* v_ptr;\
        if ( (v_ptr = object_cast<std::vector<typename AK::X> >(& (lo.approx()))) ) { \
          std::vector<typename LK::X> V;\
          V.resize(v_ptr->size());                           \
          for (unsigned int i = 0; i < v_ptr->size(); i++) {               \
            V[i] = typename LK::X(new Lazy_rep_1<typename AK::X, typename EK::X, Ith_for_intersection<typename AK::X>, \
                                                 Ith_for_intersection<typename EK::X>, E2A, Lazy_object> \
                                  (Ith_for_intersection<typename AK::X>(i), Ith_for_intersection<typename EK::X>(i), lo)); \
          }                                                           \
          return make_object(V);                                      \
        }\
      }

CGAL_Kernel_obj(Point_2)
CGAL_Kernel_obj(Point_3)  
#undef CGAL_Kernel_obj

      std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's (#1)" << std::endl;
      std::cerr << "dynamic type of the Object : " << lo.approx().type().name() << std::endl;

    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      ET eto = ec(CGAL::exact(l1), CGAL::exact(l2));
      return make_lazy<LK>(eto);
    }
    return Object();
  }

  template <typename L1, typename L2, typename L3>
  result_type
  operator()(const L1& l1, const L2& l2, const L3& l3) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      Lazy_object lo(new Lazy_rep_3<result_type, result_type, AC, EC, E2A, L1, L2, L3>(ac, ec, l1, l2, l3));

      if(lo.approx().is_empty())
        return Object();

#define CGAL_Kernel_obj(X) \
      if (object_cast<typename AK::X>(& (lo.approx()))) { \
	typedef Lazy_rep_1<typename AK::X, typename EK::X, Object_cast<typename AK::X>, Object_cast<typename EK::X>, E2A, Lazy_object> Lcr; \
	Lcr * lcr = new Lcr(Object_cast<typename AK::X>(), Object_cast<typename EK::X>(), lo); \
	return make_object(typename LK::X(lcr)); \
      }

#include <CGAL/Kernel/interface_macros.h>

      std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's (#1)" << std::endl;
      std::cerr << "dynamic type of the Object : " << lo.approx().type().name() << std::endl;

    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);
      ET eto = ec(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3));
      return make_lazy<LK>(eto);
    }
    return Object();
  }

};



//____________________________________________________________
// The magic functor that has Lazy<Something> as result type.
// Two versions are distinguished: one that needs to fiddle 
// with result_of and another that can forward the result types.

namespace internal {
BOOST_MPL_HAS_XXX_TRAIT_DEF(result_type)

// lift boost::get into a functor with a result_type member name and
// extend it to operate on optionals

// TODO there is a mismatch between the result_type typedef and the
// actual return type of operator()
template<typename T>
struct Variant_cast {
  typedef T result_type;

  template<BOOST_VARIANT_ENUM_PARAMS(typename U)>
  const T&
  operator()(const boost::optional< boost::variant< BOOST_VARIANT_ENUM_PARAMS(U) > >& o) const {
    // can throw but should never because we always build it inside
    // a static visitor with the right type
    return boost::get<T>(*o);
  }

  template<BOOST_VARIANT_ENUM_PARAMS(typename U)>
  T&
  operator()(boost::optional< boost::variant< BOOST_VARIANT_ENUM_PARAMS(U) > >& o) const {
    // can throw but should never because we always build it inside
    // a static visitor with the right type, if it throws bad_get 
    return boost::get<T>(*o);
  }
};


template<typename Result, typename AK, typename LK, typename EK, typename Origin>
struct Fill_lazy_variant_visitor_2 : boost::static_visitor<> {
  Fill_lazy_variant_visitor_2(Result& r, Origin& o) : r(&r), o(&o) {}
  Result* r;
  Origin* o;
    
  template<typename T>
  void operator()(const T&) {
    // the equivalent type we are currently matching in the lazy kernel
    typedef T AKT;
    typedef typename Type_mapper<AKT, AK, EK>::type EKT;
    typedef typename Type_mapper<AKT, AK, LK>::type LKT;

    typedef Lazy_rep_1<AKT, EKT, Variant_cast<AKT>, Variant_cast<EKT>, typename LK::E2A, Origin> Lcr;
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
      V[i] = LKT(new Lazy_rep_1<AKT, EKT, Ith_for_intersection<AKT>,
                 Ith_for_intersection<EKT>, typename LK::E2A, Origin>
                 (Ith_for_intersection<AKT>(i), Ith_for_intersection<EKT>(i), *o));
    }
      
    *r = V;
  }
};

template<typename Result, typename AK, typename LK, typename EK>
struct Fill_lazy_variant_visitor_0 : boost::static_visitor<> {
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

template <typename LK, typename AC, typename EC>
struct Lazy_construction_variant {
  static const bool Protection = true;

  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename EK::FT EFT;
  typedef typename LK::E2A E2A;


  template<typename>
  struct result {
    // this does not default, if you want to make a lazy lazy-kernel,
    // you are on your own
  };

  #define CGAL_RESULT(z, n, d) \
    template< typename F, BOOST_PP_ENUM_PARAMS(n, class T) >            \
    struct result<F( BOOST_PP_ENUM_PARAMS(n, T) )> {                    \
      BOOST_PP_REPEAT(n, CGAL_TYPEMAP_AC, T)                            \
      typedef typename Type_mapper<                                     \
        typename cpp11::result_of<AC( BOOST_PP_ENUM_PARAMS(n, A) )>::type, AK, LK>::type type; \
    };

  BOOST_PP_REPEAT_FROM_TO(1, 9, CGAL_RESULT, _)

  template <typename L1, typename L2>
  typename result<Lazy_construction_variant(L1, L2)>::type
  operator()(const L1& l1, const L2& l2) const {
    typedef typename cpp11::result_of<Lazy_construction_variant(L1, L2)>::type result_type;
    
    typedef typename cpp11::result_of<AC(typename Type_mapper<L1, LK, AK>::type, 
                                         typename Type_mapper<L2, LK, AK>::type)>::type AT;
    typedef typename cpp11::result_of<EC(typename Type_mapper<L1, LK, EK>::type, 
                                         typename Type_mapper<L2, LK, EK>::type)>::type ET;
    
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;

    try {
      Lazy<AT, ET, EFT, E2A> lazy(new Lazy_rep_2<AT, ET, AC, EC, E2A, L1, L2>(AC(), EC(), l1, l2));

      // the approximate result requires the trait with types from the AK 
      AT approx_v = lazy.approx();
      // the result we build
      result_type res;

      if(!approx_v) {
        // empty
        return res;
      }

      // the static visitor fills the result_type with the correct unwrapped type
      internal::Fill_lazy_variant_visitor_2< result_type, AK, LK, EK, Lazy<AT, ET, EFT, E2A> > visitor(res, lazy);
      boost::apply_visitor(visitor, *approx_v);
      
      return res;
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);

      ET exact_v = EC()(CGAL::exact(l1), CGAL::exact(l2));
      result_type res;

      if(!exact_v) {
        return res;
      }

      internal::Fill_lazy_variant_visitor_0<result_type, AK, LK, EK> visitor(res);
      boost::apply_visitor(visitor, *exact_v);
      return res;
    }
  }

  template <typename L1, typename L2, typename L3>
  typename result<Lazy_construction_variant(L1, L2, L3)>::type
  operator()(const L1& l1, const L2& l2, const L3& l3) const {
    typedef typename result<Lazy_construction_variant(L1, L2, L3)>::type result_type;
    
    typedef typename cpp11::result_of<AC(typename Type_mapper<L1, LK, AK>::type, 
                                         typename Type_mapper<L2, LK, AK>::type,
                                         typename Type_mapper<L3, LK, AK>::type)>::type AT;
    typedef typename cpp11::result_of<EC(typename Type_mapper<L1, LK, EK>::type, 
                                         typename Type_mapper<L2, LK, EK>::type,
                                         typename Type_mapper<L3, LK, EK>::type)>::type ET;

    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    Protect_FPU_rounding<Protection> P;
    try {
      Lazy<AT, ET, EFT, E2A> lazy(new Lazy_rep_3<AT, ET, AC, EC, E2A, L1, L2, L3>(AC(), EC(), l1, l2, l3));

      // the approximate result requires the trait with types from the AK 
      AT approx_v = lazy.approx();
      // the result we build
      result_type res;

      if(!approx_v) {
        // empty
        return res;
      }

      // the static visitor fills the result_type with the correct unwrapped type
      internal::Fill_lazy_variant_visitor_2< result_type, AK, LK, EK, Lazy<AT, ET, EFT, E2A> > visitor(res, lazy);
      boost::apply_visitor(visitor, *approx_v);
      
      return res;
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);

      ET exact_v = EC()(CGAL::exact(l1), CGAL::exact(l2), CGAL::exact(l3));
      result_type res;

      if(!exact_v) {
        return res;
      }

      internal::Fill_lazy_variant_visitor_0< result_type, AK, LK, EK> visitor(res);
      boost::apply_visitor(visitor, *exact_v);
      return res;
    }
  }
};

template<typename LK, typename AC, typename EC, typename E2A = Default, 
         bool has_result_type = internal::has_result_type<AC>::value && internal::has_result_type<EC>::value >
struct Lazy_construction;


// we have a result type, low effort
template<typename LK, typename AC, typename EC, typename E2A_>
struct Lazy_construction<LK, AC, EC, E2A_, true> {
  static const bool Protection = true;

  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename boost::remove_cv< 
    typename boost::remove_reference < typename AC::result_type >::type >::type AT;
  typedef typename boost::remove_cv< 
    typename boost::remove_reference < typename EC::result_type >::type >::type  ET;

  typedef typename EK::FT EFT;
  typedef typename Default::Get<E2A_, typename LK::E2A>::type E2A;
  
  typedef typename Type_mapper<AT, AK, LK>::type result_type;

  AC ac;
  EC ec;

#define CGAL_CONSTRUCTION_OPERATOR(z, n, d  )                           \
  template<BOOST_PP_ENUM_PARAMS(n, class L)>                            \
  result_type                                                           \
  operator()( BOOST_PP_ENUM(n, CGAL_LARGS, _) ) const {                 \
    typedef Lazy< AT, ET, EFT, E2A> Handle;                             \
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp); \
    Protect_FPU_rounding<Protection> P;                                 \
    try {                                                               \
      return result_type( Handle(new Lazy_rep_##n<AT, ET, AC, EC, E2A, BOOST_PP_ENUM_PARAMS(n, L)>(ac, ec, BOOST_PP_ENUM_PARAMS(n, l)))); \
    } catch (Uncertain_conversion_exception) {                          \
      CGAL_BRANCH_PROFILER_BRANCH(tmp);                                 \
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);          \
      return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>(ec( BOOST_PP_ENUM(n, CGAL_LEXACT, _) ))) ); \
    }                                                                   \
  }        

  // arity 1-8 
  BOOST_PP_REPEAT_FROM_TO(1, 9, CGAL_CONSTRUCTION_OPERATOR, _)

  // nullary
  result_type
  operator()() const
  {
    typedef Lazy<AT, ET, EFT, E2A> Handle;
    return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>()) );
  }

#undef CGAL_CONSTRUCTION_OPERATOR
  
};


template <typename LK, typename AC, typename EC, typename E2A_>
struct Lazy_construction<LK, AC, EC, E2A_, false>
{
  static const bool Protection = true;

  typedef typename LK::Approximate_kernel AK;
  typedef typename LK::Exact_kernel EK;
  typedef typename EK::FT EFT;
  typedef typename Default::Get<E2A_, typename LK::E2A>::type E2A;

  template<typename>
  struct result {
    // this does not default, if you want to make a lazy lazy-kernel,
    // you are on your own
  };

  AC ac;
  EC ec;

  // acquire the result_type of the approximate kernel, map it back to the lazy kernel object
#define CGAL_RESULT(z, n, d) \
template< typename F, BOOST_PP_ENUM_PARAMS(n, class T) > \
struct result<F( BOOST_PP_ENUM_PARAMS(n, T) )> { \
  BOOST_PP_REPEAT(n, CGAL_TYPEMAP_AC, T)                                     \
  typedef typename Type_mapper< typename cpp11::result_of<AC( BOOST_PP_ENUM_PARAMS(n, A) )>::type, AK, LK>::type type; \
};

  BOOST_PP_REPEAT_FROM_TO(1, 9, CGAL_RESULT, _)

#define CGAL_CONSTRUCTION_OPERATOR(z, n, d)                                      \
  template<BOOST_PP_ENUM_PARAMS(n, class L)>                            \
  typename cpp11::result_of<Lazy_construction(BOOST_PP_ENUM_PARAMS(n, L))>::type \
  operator()( BOOST_PP_ENUM(n, CGAL_LARGS, _) ) {                            \
    BOOST_PP_REPEAT(n, CGAL_TYPEMAP_EC, L)                                     \
    BOOST_PP_REPEAT(n, CGAL_TYPEMAP_AC, L)                                     \
    typedef typename boost::remove_cv< typename boost::remove_reference < \
                                        typename cpp11::result_of< EC(BOOST_PP_ENUM_PARAMS(n, E)) >::type >::type >::type ET; \
    typedef typename boost::remove_cv< typename boost::remove_reference < \
                                        typename cpp11::result_of< AC(BOOST_PP_ENUM_PARAMS(n, A)) >::type >::type >::type AT; \
    typedef Lazy< AT, ET, EFT, E2A> Handle; \
    typedef typename cpp11::result_of<Lazy_construction(BOOST_PP_ENUM_PARAMS(n, L))>::type result_type; \
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp); \
    Protect_FPU_rounding<Protection> P;                                   \
    try {                                                                 \
      return result_type( Handle(new Lazy_rep_##n<AT, ET, AC, EC, E2A, BOOST_PP_ENUM_PARAMS(n, L)>(ac, ec, BOOST_PP_ENUM_PARAMS(n, l)))); \
    } catch (Uncertain_conversion_exception) {                          \
      CGAL_BRANCH_PROFILER_BRANCH(tmp);                                 \
      Protect_FPU_rounding<!Protection> P2(CGAL_FE_TONEAREST);          \
      return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>(ec( BOOST_PP_ENUM(n, CGAL_LEXACT, _) ))) ); \
    }                                                                   \
  }        

  // arity 1-8 
  BOOST_PP_REPEAT_FROM_TO(1, 9, CGAL_CONSTRUCTION_OPERATOR, _)

  // nullary
  typename Type_mapper< typename cpp11::result_of<AC()>::type ,AK, LK>::type
  operator()() const
  {
    typedef typename cpp11::result_of<AC()>::type AT;
    typedef typename cpp11::result_of<EC()>::type ET;
    typedef Lazy<AT, ET, EFT, E2A> Handle;
    typedef typename Type_mapper< typename cpp11::result_of<AC()>::type ,AK, LK>::type result_type;

    return result_type( Handle(new Lazy_rep_0<AT,ET,E2A>()) );
  }
};

} //namespace CGAL

#undef CGAL_TYPEMAP_AC
#undef CGAL_TYPEMAP_EC
#undef CGAL_LEXACT
#undef CGAL_LARGS


#endif // CGAL_LAZY_H
