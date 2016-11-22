// Copyright (c) 1999-2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_LAZY_EXACT_NT_H
#define CGAL_LAZY_EXACT_NT_H

#define CGAL_int(T)    typename First_if_different<int,    T>::Type
#define CGAL_double(T) typename First_if_different<double, T>::Type
#define CGAL_To_interval(T) To_interval<T>


#include <CGAL/number_type_basic.h>
#include <CGAL/assertions.h>

#include <boost/iterator/transform_iterator.hpp> // for Root_of functor
#include <boost/operators.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>

#include <boost/mpl/if.hpp>
#include <boost/mpl/logical.hpp>

#include <CGAL/Interval_nt.h>
#include <CGAL/Handle.h>
#include <CGAL/NT_converter.h>

#include <CGAL/Profile_counter.h>
#include <CGAL/Lazy.h>

#include <CGAL/Sqrt_extension_fwd.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/tss.h>

#include <CGAL/IO/io.h>


/*
 * This file contains the definition of the number type Lazy_exact_nt<ET>,
 * where ET is an exact number type (must provide the exact operations needed).
 *
 * Lazy_exact_nt<ET> provides a DAG-based lazy evaluation, like LEDA's real,
 * Core's Expr, LEA's lazy rationals...
 *
 * The values are first approximated using Interval_base.
 * The exactness is provided when needed by ET.
 *
 * Lazy_exact_nt<ET> is just a handle to the abstract base class
 * Lazy_exact_nt_rep which has pure virtual methods .approx() and .exact().
 * From this class derives one class per operation, with one constructor.
 *
 * The DAG is managed by :
 * - Handle and Rep.
 * - virtual functions to denote the various operators (instead of an enum).
 *
 * Other packages with vaguely similar design : APU, MetaCGAL, LOOK.
 */

/*
 * TODO :
 * - Generalize it for constructions at the kernel level.
 * - Add mixed operations with ET too ?
 * - Interval refinement functionnality ?
 * - Separate the handle and the representation(s) in 2 files (?)
 *   maybe not a good idea, better if everything related to one operation is
 *   close together.
 * - Add a CT template parameter ?
 * - Add a string constant to provide an expression string (a la MetaCGAL) ?
 *   // virtual ostream operator<<() const = 0; // or string, like Core ?
 * - Have a template-expression (?) thing that evaluates a temporary element,
 *   and allocates stuff in memory only when really needs to convert to the
 *   NT.  (cf gmp++, and maybe other things, Blitz++, Synaps...)
 */

/*
 * Interface of the rep classes:
 * - .approx()      returns Interval_nt<> (assumes rounding=nearest).
 *                  [ only called from the handle, and declared in the base ]
 * - .exact()       returns ET, if not already done, computes recursively
 *
 * - .rafine_approx()   ??
 */

namespace CGAL {

template <class NT> class Lazy_exact_nt;


#ifdef CGAL_LAZY_KERNEL_DEBUG
template <typename ET>
inline
void
print_dag(const Lazy_exact_nt<ET>& l, std::ostream& os, int level=0)
{
  l.print_dag(os, level);
}
#endif

// Abstract base representation class for lazy numbers
template <typename ET>
struct Lazy_exact_nt_rep : public Lazy_exact_nt<ET>::Self_rep
{
  typedef typename Lazy_exact_nt<ET>::Self_rep  Base;

  Lazy_exact_nt_rep (const Interval_nt<false> & i)
      : Base(i) {}

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
  }
#endif
};

// int constant
template <typename ET>
struct Lazy_exact_Int_Cst : public Lazy_exact_nt_rep<ET>
{
  Lazy_exact_Int_Cst (int i)
      : Lazy_exact_nt_rep<ET>(double(i)) {}

  void update_exact() const { this->et = new ET((int)this->approx().inf()); }
};

// double constant
template <typename ET, typename X>
struct Lazy_exact_Cst : public Lazy_exact_nt_rep<ET>
{
  Lazy_exact_Cst (X x)
      : Lazy_exact_nt_rep<ET>(x), cste(x) {}

  void update_exact() const { this->et = new ET(cste); }

  private:
  X cste;
};

// Exact constant
template <typename ET>
struct Lazy_exact_Ex_Cst : public Lazy_exact_nt_rep<ET>
{
  Lazy_exact_Ex_Cst (const ET & e)
      : Lazy_exact_nt_rep<ET>(CGAL_NTS to_interval(e))
  {
    this->et = new ET(e);
  }

  void update_exact() const { CGAL_error(); }
};

// Construction from a Lazy_exact_nt<ET1> (which keeps the lazyness).
template <typename ET, typename ET1>
class Lazy_lazy_exact_Cst : public Lazy_exact_nt_rep<ET>
{
  mutable Lazy_exact_nt<ET1> l;

public:

  Lazy_lazy_exact_Cst (const Lazy_exact_nt<ET1> &x)
      : Lazy_exact_nt_rep<ET>(x.approx()), l(x)
  {
    this->set_depth(l.depth() + 1);
  }

  void update_exact() const
  {
    this->et = new ET(l.exact());
    this->at = l.approx();
    prune_dag();
  }

  void prune_dag() const { l = Lazy_exact_nt<ET1>(); }
};


// Unary  operations: abs, sqrt, square.
// Binary operations: +, -, *, /, min, max.

// Base unary operation
template <typename ET>
struct Lazy_exact_unary : public Lazy_exact_nt_rep<ET>
{
  mutable Lazy_exact_nt<ET> op1;

  Lazy_exact_unary (const Interval_nt<false> &i, const Lazy_exact_nt<ET> &a)
      : Lazy_exact_nt_rep<ET>(i), op1(a)
  {
    this->set_depth(op1.depth() + 1);
  }

  void prune_dag() const { op1 = Lazy_exact_nt<ET>(); }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      msg(os, level, "Unary number operator:");
      CGAL::print_dag(op1, os, level+1);
    }
  }
#endif
};

// Base binary operation
template <typename ET, typename ET1 = ET, typename ET2 = ET>
struct Lazy_exact_binary : public Lazy_exact_nt_rep<ET>
{
  mutable Lazy_exact_nt<ET1> op1;
  mutable Lazy_exact_nt<ET2> op2;

  Lazy_exact_binary (const Interval_nt<false> &i,
		     const Lazy_exact_nt<ET1> &a, const Lazy_exact_nt<ET2> &b)
      : Lazy_exact_nt_rep<ET>(i), op1(a), op2(b)
  {
    this->set_depth((std::max)(op1.depth(), op2.depth()) + 1);
  }

  void prune_dag() const
  {
    op1 = Lazy_exact_nt<ET1>();
    op2 = Lazy_exact_nt<ET2>();
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      msg(os, level, "Binary number operator:");
      CGAL::print_dag(op1, os, level+1);
      CGAL::print_dag(op2, os, level+1);
    }
  }
#endif
};

// Here we could use a template class for all operations (STL provides
// function objects plus, minus, multiplies, divides...).  But it would require
// a template template parameter, and GCC 2.95 seems to crash easily with them.

// Macro for unary operations
#define CGAL_LAZY_UNARY_OP(OP, NAME)                                     \
template <typename ET>                                                   \
struct NAME : public Lazy_exact_unary<ET>                                \
{                                                                        \
  typedef typename Lazy_exact_unary<ET>::AT::Protector P;                \
  NAME (const Lazy_exact_nt<ET> &a)                                      \
      : Lazy_exact_unary<ET>((P(), OP(a.approx())), a) {}                \
                                                                         \
  void update_exact() const                                              \
  {                                                                      \
    this->et = new ET(OP(this->op1.exact()));                            \
    if (!this->approx().is_point())                                      \
      this->at = CGAL_NTS to_interval(*(this->et));                      \
    this->prune_dag();                                                   \
   }                                                                     \
};

CGAL_LAZY_UNARY_OP(opposite,  Lazy_exact_Opp)
CGAL_LAZY_UNARY_OP(CGAL_NTS abs,    Lazy_exact_Abs)
CGAL_LAZY_UNARY_OP(CGAL_NTS square, Lazy_exact_Square)
CGAL_LAZY_UNARY_OP(CGAL_NTS sqrt,   Lazy_exact_Sqrt)

// A macro for +, -, * and /
#define CGAL_LAZY_BINARY_OP(OP, NAME)                                    \
template <typename ET, typename ET1 = ET, typename ET2 = ET>             \
struct NAME : public Lazy_exact_binary<ET, ET1, ET2>                     \
{                                                                        \
  typedef typename Lazy_exact_binary<ET,ET1,ET2>::AT::Protector P;	 \
  NAME (const Lazy_exact_nt<ET1> &a, const Lazy_exact_nt<ET2> &b)        \
    : Lazy_exact_binary<ET, ET1, ET2>((P(), a.approx() OP b.approx()), a, b) {} \
                                                                         \
  void update_exact() const                                              \
  {                                                                      \
    this->et = new ET(this->op1.exact() OP this->op2.exact());           \
    if (!this->approx().is_point())                                      \
      this->at = CGAL_NTS to_interval(*(this->et));                      \
    this->prune_dag();                                                   \
   }                                                                     \
};

CGAL_LAZY_BINARY_OP(+, Lazy_exact_Add)
CGAL_LAZY_BINARY_OP(-, Lazy_exact_Sub)
CGAL_LAZY_BINARY_OP(*, Lazy_exact_Mul)
CGAL_LAZY_BINARY_OP(/, Lazy_exact_Div)

// Minimum
template <typename ET>
struct Lazy_exact_Min : public Lazy_exact_binary<ET>
{
  Lazy_exact_Min (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_binary<ET>((CGAL::min)(a.approx(), b.approx()), a, b) {}

  void update_exact() const
  {
    this->et = new ET((CGAL::min)(this->op1.exact(), this->op2.exact()));
    if (!this->approx().is_point()) 
      this->at = CGAL_NTS to_interval(*(this->et));
    this->prune_dag();
  }
};

// Maximum
template <typename ET>
struct Lazy_exact_Max : public Lazy_exact_binary<ET>
{
  Lazy_exact_Max (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_binary<ET>((CGAL::max)(a.approx(), b.approx()), a, b) {}

  void update_exact() const
  {
    this->et = new ET((CGAL::max)(this->op1.exact(), this->op2.exact()));
    if (!this->approx().is_point()) 
      this->at = CGAL_NTS to_interval(*(this->et));
    this->prune_dag();
  }
};


// The real number type, handle class
template <typename ET_>
class Lazy_exact_nt
  : public Lazy<Interval_nt<false>, ET_, Lazy_exact_nt<ET_>, To_interval<ET_> >
  , boost::ordered_euclidian_ring_operators2< Lazy_exact_nt<ET_>, int >
  , boost::ordered_euclidian_ring_operators2< Lazy_exact_nt<ET_>, double >
{
public:

  typedef Lazy_exact_nt<ET_> Self;
  typedef Lazy<Interval_nt<false>, ET_, Self, To_interval<ET_> > Base;
  typedef typename Base::Self_rep  Self_rep;

  typedef typename Base::ET ET; // undocumented
  typedef typename Base::AT AT; // undocumented

  typedef typename Base::Exact_type        Exact_type;
  typedef typename Base::Approximate_type  Approximate_type;

public :

  Lazy_exact_nt () {}

  Lazy_exact_nt (Self_rep *r)
    : Base(r) {}

  // Also check that ET and AT are constructible from T?
  template<class T>
  Lazy_exact_nt (T i, typename boost::enable_if<boost::mpl::and_<
      boost::mpl::or_<boost::is_arithmetic<T>, boost::is_enum<T> >,
      boost::mpl::not_<boost::is_same<T,ET> > >,void*>::type=0)
    : Base(new Lazy_exact_Cst<ET,T>(i)) {}

  Lazy_exact_nt (const ET & e)
    : Base(new Lazy_exact_Ex_Cst<ET>(e)){}

  template <class ET1>
  Lazy_exact_nt (const Lazy_exact_nt<ET1> &x,
      typename boost::enable_if<is_implicit_convertible<ET1,ET>,int>::type=0)
    : Base(new Lazy_lazy_exact_Cst<ET, ET1>(x)){}

  template <class ET1>
  explicit Lazy_exact_nt (const Lazy_exact_nt<ET1> &x,
  typename boost::disable_if<is_implicit_convertible<ET1,ET>,int>::type=0)
    : Base(new Lazy_lazy_exact_Cst<ET, ET1>(x)){}

  Self operator+ () const
  { return *this; }

  Self operator- () const
  { return new Lazy_exact_Opp<ET>(*this); }

  Self & operator+=(const Self& b)
  { return *this = new Lazy_exact_Add<ET>(*this, b); }

  Self & operator-=(const Self& b)
  { return *this = new Lazy_exact_Sub<ET>(*this, b); }

  Self & operator*=(const Self& b)
  { return *this = new Lazy_exact_Mul<ET>(*this, b); }

  Self & operator/=(const Self& b)
  {
    CGAL_precondition(b != 0);
    return *this = new Lazy_exact_Div<ET>(*this, b);
  }

  // Mixed operators. (could be optimized ?)
  Self & operator+=(CGAL_int(ET) b)
  { return *this = new Lazy_exact_Add<ET>(*this, b); }

  Self & operator-=(CGAL_int(ET) b)
  { return *this = new Lazy_exact_Sub<ET>(*this, b); }

  Self & operator*=(CGAL_int(ET) b)
  { return *this = new Lazy_exact_Mul<ET>(*this, b); }

  Self & operator/=(CGAL_int(ET) b)
  {
    CGAL_precondition(b != 0);
    return *this = new Lazy_exact_Div<ET>(*this, b);
  }

  Self & operator+=(CGAL_double(ET) b)
  { return *this = new Lazy_exact_Add<ET>(*this, b); }

  Self & operator-=(CGAL_double(ET) b)
  { return *this = new Lazy_exact_Sub<ET>(*this, b); }

  Self & operator*=(CGAL_double(ET) b)
  { return *this = new Lazy_exact_Mul<ET>(*this, b); }

  Self & operator/=(CGAL_double(ET) b)
  {
    CGAL_precondition(b != 0);
    return *this = new Lazy_exact_Div<ET>(*this, b);
  }

  // % kills filtering
  Self & operator%=(const Self& b)
  {
    CGAL_precondition(b != 0);
    ET res = this->exact();
    res %= b.exact();
    return *this = Lazy_exact_nt<ET>(res);
  }

  Self & operator%=(int b)
  {
    CGAL_precondition(b != 0);
    ET res = this->exact();
    res %= b;
    return *this = Lazy_exact_nt<ET>(res);
  }

  Interval_nt<true> interval() const
  {
    const Interval_nt<false>& i = this->approx();
    return Interval_nt<true>(i.inf(), i.sup());
  }

  Interval_nt_advanced approx_adv() const
  { return this->ptr()->approx(); }

private:
  static double & relative_precision_of_to_double_internal()
  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(double, relative_precision_of_to_double, 0.00001);
      return relative_precision_of_to_double;
  }

public:
  static const double & get_relative_precision_of_to_double()
  {
    return relative_precision_of_to_double_internal();
  }

  static void set_relative_precision_of_to_double(const double & d)
  {
      CGAL_assertion((0 < d) & (d < 1));
      relative_precision_of_to_double_internal() = d;
  }

  bool identical(const Self& b) const
  {
      return ::CGAL::identical(
              static_cast<const Handle &>(*this),
              static_cast<const Handle &>(b));
  }

  template < typename T >
  bool identical(const T&) const
  { return false; }

};



template <typename ET1, typename ET2>
bool
operator<(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{
  CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
  if (a.identical(b))
    return false;
  Uncertain<bool> res = a.approx() < b.approx();
  if (is_certain(res))
    return get_certain(res);
  CGAL_BRANCH_PROFILER_BRANCH(tmp);
  return a.exact() < b.exact();
}

template <typename ET1, typename ET2>
bool
operator==(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{
  CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
  if (a.identical(b))
    return true;
  Uncertain<bool> res = a.approx() == b.approx();
  if (is_certain(res))
    return get_certain(res);
  CGAL_BRANCH_PROFILER_BRANCH(tmp);
  return a.exact() == b.exact();
}

template <typename ET1, typename ET2>
inline
bool
operator>(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{ return b < a; }

template <typename ET1, typename ET2>
inline
bool
operator>=(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{ return ! (a < b); }

template <typename ET1, typename ET2>
inline
bool
operator<=(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{ return b >= a; }

template <typename ET1, typename ET2>
inline
bool
operator!=(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{ return ! (a == b); }


template <typename ET>
inline
Lazy_exact_nt<ET>
operator%(const Lazy_exact_nt<ET>& a, const Lazy_exact_nt<ET>& b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  CGAL_precondition(b != 0);
  return Lazy_exact_nt<ET>(a) %= b;
}



// Mixed operators with int.
template <typename ET>
bool
operator<(const Lazy_exact_nt<ET>& a, int b)
{
  CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
  Uncertain<bool> res = a.approx() < b;
  if (is_certain(res))
    return res;
  CGAL_BRANCH_PROFILER_BRANCH(tmp);
  return a.exact() < b;
}

template <typename ET>
bool
operator>(const Lazy_exact_nt<ET>& a, int b)
{
  CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
  Uncertain<bool> res = b < a.approx();
  if (is_certain(res))
    return get_certain(res);
  CGAL_BRANCH_PROFILER_BRANCH(tmp);
  return b < a.exact();
}

template <typename ET>
bool
operator==(const Lazy_exact_nt<ET>& a, int b)
{
  CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
  Uncertain<bool> res = b == a.approx();
  if (is_certain(res))
    return get_certain(res);
  CGAL_BRANCH_PROFILER_BRANCH(tmp);
  return b == a.exact();
}


// Mixed operators with double.
template <typename ET>
bool
operator<(const Lazy_exact_nt<ET>& a, double b)
{
  CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
  Uncertain<bool> res = a.approx() < b;
  if (is_certain(res))
    return res;
  CGAL_BRANCH_PROFILER_BRANCH(tmp);
  return a.exact() < b;
}

template <typename ET>
bool
operator>(const Lazy_exact_nt<ET>& a, double b)
{
  CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
  Uncertain<bool> res = b < a.approx();
  if (is_certain(res))
    return res;
  CGAL_BRANCH_PROFILER_BRANCH(tmp);
  return b < a.exact();
}

template <typename ET>
bool
operator==(const Lazy_exact_nt<ET>& a, double b)
{
  CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
  Uncertain<bool> res = b == a.approx();
  if (is_certain(res))
    return res;
  CGAL_BRANCH_PROFILER_BRANCH(tmp);
  return b == a.exact();
}



template <typename ET1, typename ET2>
Lazy_exact_nt< typename Coercion_traits<ET1, ET2>::Type >
operator+(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  return new Lazy_exact_Add<typename Coercion_traits<ET1, ET2>::Type,
                            ET1, ET2>(a, b);
}

template <typename ET1, typename ET2>
Lazy_exact_nt< typename Coercion_traits<ET1, ET2>::Type >
operator-(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  return new Lazy_exact_Sub<typename Coercion_traits<ET1, ET2>::Type,
                            ET1, ET2>(a, b);
}

template <typename ET1, typename ET2>
Lazy_exact_nt< typename Coercion_traits<ET1, ET2>::Type >
operator*(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  return new Lazy_exact_Mul<typename Coercion_traits<ET1, ET2>::Type,
                            ET1, ET2>(a, b);
}

template <typename ET1, typename ET2>
Lazy_exact_nt< typename Coercion_traits<ET1, ET2>::Type >
operator/(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  CGAL_precondition(b != 0);
  return new Lazy_exact_Div<typename Coercion_traits<ET1, ET2>::Type,
                            ET1, ET2>(a, b);
}

//
// Algebraic structure traits
//

namespace INTERN_LAZY_EXACT_NT {

template< class NT, class Functor >
struct Simplify_selector {
  struct Simplify : public std::unary_function<NT&, void> {
    void operator()( NT& ) const {
      // TODO: In the old implementation the Simplify-functor was the default
      //       (which does nothing). But this cannot be the correct way!?
    }
  };
};

template< class NT >
struct Simplify_selector< NT, Null_functor > {
  typedef Null_functor Simplify;
};

template< class NT, class Functor >
struct Unit_part_selector {
  struct Unit_part : public std::unary_function<NT, NT > {
    NT operator()( const NT& x ) const {
      return NT( CGAL_NTS unit_part( x.exact() ) );
    }
  };
};

template< class NT >
struct Unit_part_selector< NT, Null_functor > {
  typedef Null_functor Unit_part;
};

template< class NT, class Functor >
struct Is_zero_selector {
  struct Is_zero : public std::unary_function<NT, bool > {
    bool operator()( const NT& x ) const {
      return CGAL_NTS is_zero( x.exact() );
    }
  };
};

template< class NT >
struct Is_zero_selector< NT, Null_functor > {
  typedef Null_functor Is_zero;
};

template< class NT, class Functor >
struct Is_one_selector {
  struct Is_one : public std::unary_function<NT, bool > {
    bool operator()( const NT& x ) const {
      return CGAL_NTS is_one( x.exact() );
    }
  };
};

template< class NT >
struct Is_one_selector< NT, Null_functor > {
  typedef Null_functor Is_one;
};

template< class NT, class Functor >
struct Square_selector {
  struct Square : public std::unary_function<NT, NT > {
    NT operator()( const NT& x ) const {
      CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
      return new Lazy_exact_Square<typename NT::ET>(x);
    }
  };
};

template< class NT >
struct Square_selector< NT, Null_functor > {
  typedef Null_functor Square;
};

template< class NT, class Functor >
struct Integral_division_selector {
  struct Integral_division : public std::binary_function<NT, NT, NT > {
    NT operator()( const NT& x, const NT& y ) const {
      return NT( CGAL_NTS integral_division( x.exact(), y.exact() ) );
    }

    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( NT )
  };
};

template< class NT >
struct Integral_division_selector< NT, Null_functor > {
  typedef Null_functor Integral_division;
};

template< class NT, class Functor >
struct Is_square_selector {
  struct Is_square : public std::binary_function<NT, NT&, bool > {
      bool operator()( const NT& x, NT& y ) const {
          typename NT::ET y_et;
          bool result = CGAL_NTS is_square( x.exact(), y_et );
          y = NT(y_et);
          return result;
      }
      bool operator()( const NT& x) const {
          typename NT::ET y_et;
          return CGAL_NTS is_square( x.exact(), y_et );
      }
  };
};

template< class NT >
struct Is_square_selector< NT, Null_functor > {
  typedef Null_functor Is_square;
};


template <class NT, class AlgebraicStructureTag>
struct Sqrt_selector{
    struct Sqrt : public std::unary_function<NT, NT > {
        NT operator ()(const NT& x) const {
          CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
          CGAL_precondition(x >= 0);
          return new Lazy_exact_Sqrt<typename NT::ET>(x);
        }
    };
};
template <class NT>
struct Sqrt_selector<NT,Null_functor> {
    typedef Null_functor Sqrt;
};

template< class NT, class Functor >
struct Kth_root_selector {
  struct Kth_root : public std::binary_function<int, NT, NT > {
    NT operator()( int k, const NT& x ) const {
      return NT( CGAL_NTS kth_root( k, x.exact() ) );
    }
  };
};

template< class NT >
struct Kth_root_selector< NT, Null_functor > {
  typedef Null_functor Kth_root;
};

template< class NT, class Functor >
struct Root_of_selector {
  private:
      struct Cast{
        typedef typename NT::ET result_type;
        result_type operator()(const NT& lazy_exact) const {
          return lazy_exact.exact();
        }
      };

  public:
    struct Root_of {
//      typedef typename Functor::Boundary Boundary;
      typedef NT result_type;
      template< class Input_iterator >
      NT operator()( int k, Input_iterator begin, Input_iterator end ) const {
       Cast cast;
       return NT( typename Algebraic_structure_traits<typename NT::ET>::
                                 Root_of()( k,
                              ::boost::make_transform_iterator( begin, cast ),
                              ::boost::make_transform_iterator( end, cast ) ) );
      }

    };
};

template< class NT >
struct Root_of_selector< NT, Null_functor > {
  typedef Null_functor Root_of;
};

template< class NT, class Functor >
struct Gcd_selector {
  struct Gcd : public std::binary_function<NT, NT, NT > {
    NT operator()( const NT& x, const NT& y ) const {
     CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
     return NT( CGAL_NTS gcd( x.exact(), y.exact() ) );
    }

    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( NT )
  };
};

template< class NT >
struct Gcd_selector< NT, Null_functor > {
  typedef Null_functor Gcd;
};

template< class NT, class Functor >
struct Div_selector {
  struct Div : public std::binary_function<NT, NT, NT > {
    NT operator()( const NT& x, const NT& y ) const {
      return NT( CGAL_NTS div( x.exact(), y.exact() ) );
    }

    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( NT )

  };
};

template< class NT >
struct Div_selector< NT, Null_functor > {
  typedef Null_functor Div;
};

template< class NT, class Functor >
struct Inverse_selector {
  struct Inverse  {
    typedef NT result_type;
    NT operator()( const NT& x ) const {
      return NT( 1 ) / x ;
    }
  };
};

template< class NT >
struct Inverse_selector< NT, Null_functor > {
  typedef Null_functor Inverse;
};

template< class NT, class Functor >
struct Mod_selector {
  struct Mod : public std::binary_function<NT, NT, NT > {
    NT operator()( const NT& x, const NT& y ) const {
      return NT( CGAL_NTS mod( x.exact(), y.exact() ) );
    }

    CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( NT )

  };
};

template< class NT >
struct Mod_selector< NT, Null_functor > {
  typedef Null_functor Mod;
};

template< class NT, class Functor >
struct Div_mod_selector {
  struct Div_mod {
    typedef void result_type;
    typedef NT   first_argument_type;
    typedef NT   second_argument_type;
    typedef NT&  third_argument_type;
    typedef NT&  fourth_argument_type;

    void operator()( const NT& x, const NT& y, NT& q, NT& r ) const {
      typename NT::ET q_et;
      typename NT::ET r_et;
      CGAL_NTS div_mod( x.exact(), y.exact(), q_et, r_et );
      q = NT( q_et );
      r = NT( r_et );
    }

    template< class NT1, class NT2 >
    void operator()( const NT1& x, const NT2& y,
                     NT& q,
                     NT& r ) const {
      CGAL_static_assertion((::boost::is_same<
        typename Coercion_traits< NT1, NT2 >::Type, NT
                                              >::value));

      typename Coercion_traits< NT1, NT2 >::Cast cast;
      operator()( cast(x), cast(y), q, r );
    }
  };
};

template< class NT >
struct Div_mod_selector< NT, Null_functor >{
  typedef Null_functor Div_mod;
};
} // namespace INTERN_LAZY_EXACT_NT

template <class ET>
class Algebraic_structure_traits< Lazy_exact_nt<ET> >
    :public Algebraic_structure_traits_base
      < Lazy_exact_nt<ET>,
       typename Algebraic_structure_traits<ET>::Algebraic_category >
{
private:
    typedef Algebraic_structure_traits<ET> AST_ET;
    typedef typename AST_ET::Algebraic_category ET_as_tag;
public:
    typedef typename AST_ET::Is_exact               Is_exact;
    typedef typename AST_ET::Is_numerical_sensitive Is_numerical_sensitive;

    typedef typename INTERN_LAZY_EXACT_NT::Simplify_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Simplify > ::Simplify Simplify;

    typedef typename INTERN_LAZY_EXACT_NT::Unit_part_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Unit_part > ::Unit_part Unit_part;

    typedef typename INTERN_LAZY_EXACT_NT::Is_zero_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Is_zero > ::Is_zero Is_zero;

    typedef typename INTERN_LAZY_EXACT_NT::Is_one_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Is_one > ::Is_one Is_one;

    typedef typename INTERN_LAZY_EXACT_NT::Square_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Square > ::Square Square;

    typedef typename INTERN_LAZY_EXACT_NT::Integral_division_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Integral_division> ::Integral_division Integral_division;

    typedef typename INTERN_LAZY_EXACT_NT::Is_square_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Is_square > ::Is_square Is_square;

    typedef typename INTERN_LAZY_EXACT_NT::Sqrt_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Sqrt> ::Sqrt Sqrt;

    typedef typename INTERN_LAZY_EXACT_NT::Kth_root_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Kth_root > ::Kth_root Kth_root;

    typedef typename INTERN_LAZY_EXACT_NT::Root_of_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Root_of > ::Root_of Root_of;

    typedef typename INTERN_LAZY_EXACT_NT::Gcd_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Gcd > ::Gcd Gcd;

    typedef typename INTERN_LAZY_EXACT_NT::Div_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Div > ::Div Div;

    typedef typename INTERN_LAZY_EXACT_NT::Mod_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Mod > ::Mod Mod;

    typedef typename INTERN_LAZY_EXACT_NT::Div_mod_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Div_mod > ::Div_mod Div_mod;
    
    typedef typename INTERN_LAZY_EXACT_NT::Inverse_selector
    <Lazy_exact_nt<ET>, typename AST_ET::Inverse > ::Inverse Inverse;    
};



//
// Real embeddalbe traits
//

template < typename ET > class Real_embeddable_traits< Lazy_exact_nt<ET> >
  : public INTERN_RET::Real_embeddable_traits_base< Lazy_exact_nt<ET> , CGAL::Tag_true > {

  // Every type ET of Lazy_exact_nt<ET> has to be real embeddable.
  CGAL_static_assertion((::boost::is_same< typename Real_embeddable_traits< ET >
                                ::Is_real_embeddable, Tag_true >::value));

  public:
    typedef Lazy_exact_nt<ET> Type;

    class Abs
      : public std::unary_function< Type, Type > {
      public:
        Type operator()( const Type& a ) const {
            CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
            return new Lazy_exact_Abs<ET>(a);
        }
    };

    class Sgn
      : public std::unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& a ) const {
            CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
            Uncertain< ::CGAL::Sign> res = CGAL_NTS sign(a.approx());
            if (is_certain(res))
                return get_certain(res);
            CGAL_BRANCH_PROFILER_BRANCH(tmp);
            return CGAL_NTS sign(a.exact());
        }
    };

    class Compare
      : public std::binary_function< Type, Type,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Type& a,
                                            const Type& b ) const {
            CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
            if (a.identical(b))
                return EQUAL;
            Uncertain<Comparison_result> res = CGAL_NTS compare(a.approx(), b.approx());
            if (is_certain(res))
                return get_certain(res);
            CGAL_BRANCH_PROFILER_BRANCH(tmp);
            return CGAL_NTS compare(a.exact(), b.exact());
        }

        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Type,
                                                      Comparison_result )

    };

    class To_double
      : public std::unary_function< Type, double > {
      public:
        double operator()( const Type& a ) const {
            CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);

            const Interval_nt<false>& app = a.approx();
            double r;
            if (fit_in_double(app, r))
                return r;

            // If it's precise enough, then OK.
            if (has_smaller_relative_precision(app,
                 Lazy_exact_nt<ET>::get_relative_precision_of_to_double()))
                return CGAL_NTS to_double(app);

            CGAL_BRANCH_PROFILER_BRANCH(tmp);

            // Otherwise we trigger exact computation first,
            // which will refine the approximation.
            a.exact();
            return CGAL_NTS to_double(a.approx());
        }
    };

    class To_interval
      : public std::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& a ) const {
            CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
            return a.approx().pair();
        }
    };

    class Is_finite
      : public std::unary_function< Type, bool > {
      public:
        bool operator()( const Type& x ) const {
          return CGAL_NTS is_finite(x.approx()) || CGAL_NTS is_finite(x.exact());
        }
    };

};

template <class ET1, class ET2, class F>
class Lazy_exact_nt_coercion_traits_base {
public:
    typedef Tag_false Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
    //typedef Null_type    Type
    typedef Null_functor Cast;
};

template <class ET1, class ET2>
class Lazy_exact_nt_coercion_traits_base < Lazy_exact_nt<ET1>, Lazy_exact_nt<ET2>, Tag_true >
{
    typedef Coercion_traits<ET1,ET2> CT;
    typedef Lazy_exact_nt<ET1> A;
    typedef Lazy_exact_nt<ET2> B;
public:
    typedef Lazy_exact_nt<typename CT::Type> Type;
    typedef typename CT::Are_implicit_interoperable Are_explicit_interoperable;
    typedef typename CT::Are_implicit_interoperable Are_implicit_interoperable;

    class Cast{
    private:
        template <class T>
        Type cast(const T& x) const{ return Type(x); }
        Type cast(const Type& x) const{ return x; }
    public:
        typedef Type result_type;
        Type operator()(const A& x) const { return cast(x);}
        Type operator()(const B& x) const { return cast(x);}
    };
};


CGAL_DEFINE_COERCION_TRAITS_FOR_SELF_TEM(Lazy_exact_nt<ET>, class ET)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO_TEM(ET,Lazy_exact_nt<ET>,class ET)

template<class ET1, class ET2 >
struct Coercion_traits< Lazy_exact_nt<ET1>, Lazy_exact_nt<ET2> >
    :public Lazy_exact_nt_coercion_traits_base
           <Lazy_exact_nt<ET1>, Lazy_exact_nt<ET2>,
            typename Coercion_traits<ET1,ET2>::Are_implicit_interoperable>{};


#define CGAL_COERCION_TRAITS_LAZY_EXACT(NTX)                            \
    template<class ET>                                                  \
    struct Coercion_traits< NTX, Lazy_exact_nt<ET> >{                   \
    private:                                                            \
        typedef Coercion_traits<NTX,ET> CT;                             \
        typedef Lazy_exact_nt<ET> NT;                                   \
    public:                                                             \
        typedef typename CT::Are_explicit_interoperable                 \
        Are_explicit_interoperable;                                     \
        typedef typename CT::Are_implicit_interoperable                 \
        Are_implicit_interoperable;                                     \
    private:                                                            \
        static const  bool interoperable                                \
        =boost::is_same< Are_implicit_interoperable, Tag_false>::value; \
    public:                                                             \
        typedef typename boost::mpl::if_c <interoperable,Null_tag,NT>   \
        ::type  Type;                                          \
        typedef typename boost::mpl::if_c <interoperable, Null_functor, \
    INTERN_CT::Cast_from_to<NTX,NT> >::type Cast;                       \
    };                                                                  \
                                                                        \
    template<class ET>                                                  \
    struct Coercion_traits< Lazy_exact_nt<ET>, NTX >                    \
        :public Coercion_traits<NTX, Lazy_exact_nt<ET> >{};             \


CGAL_COERCION_TRAITS_LAZY_EXACT(int)
CGAL_COERCION_TRAITS_LAZY_EXACT(short)
CGAL_COERCION_TRAITS_LAZY_EXACT(double)
CGAL_COERCION_TRAITS_LAZY_EXACT(float)
#undef CGAL_COERCION_TRAITS_LAZY_EXACT

namespace INTERN_LAZY_EXACT_NT {

template < typename NT, typename TAG  > class Fraction_traits_base;

template < class ET >
class Fraction_traits_base <Lazy_exact_nt<ET> , CGAL::Tag_false>
    : public Fraction_traits<ET> {
public:
    typedef Lazy_exact_nt<ET>  Type;
};

template < class ET >
class Fraction_traits_base <Lazy_exact_nt<ET> , CGAL::Tag_true>{
    typedef Fraction_traits<ET> ETT;
    typedef typename ETT::Numerator_type ET_numerator;
    typedef typename ETT::Denominator_type ET_denominator;
public:
    typedef Lazy_exact_nt<ET>  Type;
    typedef Tag_true Is_fraction;
    typedef Lazy_exact_nt<ET_numerator> Numerator_type;
    typedef Lazy_exact_nt<ET_denominator> Denominator_type;

    struct Common_factor : std::binary_function<Denominator_type,Denominator_type,Denominator_type>{
        Denominator_type operator()(const Denominator_type& a, const Denominator_type& b) const {
            typename ETT::Common_factor common_factor;
            return Denominator_type(common_factor(a.exact(),b.exact()));
        }
    };
    struct Compose : std::binary_function<Type,Numerator_type,Denominator_type>{
        Type operator()(const Numerator_type& n, const Denominator_type& d) const {
            typename ETT::Compose compose;
            return Type(compose(n.exact(),d.exact()));
        }
    };
    struct Decompose {
        typedef void result_type;
        typedef Type first_argument_type;
        typedef Numerator_type second_argument_type;
        typedef Denominator_type third_argument_type;
        void operator()(const Type& f, Numerator_type& n, Denominator_type& d) const {
            typename ETT::Decompose decompose;
            ET_numerator nn;
            ET_denominator dd;
            decompose(f.exact(),nn,dd);
            n = Numerator_type(nn);
            d = Denominator_type(dd);
        }
    };
};
} // namespace INTERN_LAZY_EXACT_NT


template < class ET >
class Fraction_traits< Lazy_exact_nt< ET > >
    :public INTERN_LAZY_EXACT_NT::Fraction_traits_base<Lazy_exact_nt<ET>,
            typename Fraction_traits<ET>::Is_fraction>
{};

template < class ET >
struct Min <Lazy_exact_nt<ET> >
    : public std::binary_function<Lazy_exact_nt<ET>,Lazy_exact_nt<ET>,Lazy_exact_nt<ET> > {

    Lazy_exact_nt<ET> operator()( const Lazy_exact_nt<ET>& x, const Lazy_exact_nt<ET>& y) const
    {
      if (x.identical(y)){
        return x;
      }
      Uncertain<bool> res = x.approx() < y.approx();
      if(is_certain(res)){
        return res.make_certain() ? x : y;
      }
      CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
      return new Lazy_exact_Min<ET>(x, y);
    }
};

template < class ET >
struct Max <Lazy_exact_nt<ET> >
    : public std::binary_function<Lazy_exact_nt<ET>,Lazy_exact_nt<ET>,Lazy_exact_nt<ET> > {

    Lazy_exact_nt<ET> operator()( const Lazy_exact_nt<ET>& x, const Lazy_exact_nt<ET>& y) const
    {
      if (x.identical(y)){
        return x;
      }
      Uncertain<bool> res = x.approx() > y.approx();
      if(is_certain(res)){
        return  res.make_certain() ? x : y;
      }
        CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
        return new Lazy_exact_Max<ET>(x, y);
    }
};

template<typename ET> inline 
Lazy_exact_nt<ET> min BOOST_PREVENT_MACRO_SUBSTITUTION(
const Lazy_exact_nt<ET> & x,
const Lazy_exact_nt<ET> & y){
  return CGAL::Min<Lazy_exact_nt<ET> > ()(x,y);
}
template<typename ET> inline 
Lazy_exact_nt<ET> max BOOST_PREVENT_MACRO_SUBSTITUTION(
const Lazy_exact_nt<ET> & x,
const Lazy_exact_nt<ET> & y){
  return CGAL::Max<Lazy_exact_nt<ET> > ()(x,y);
}

template <typename ET>
std::ostream &
operator<< (std::ostream & os, const Lazy_exact_nt<ET> & a)
{ return os << CGAL_NTS to_double(a); }

template <typename ET>
std::istream &
operator>> (std::istream & is, Lazy_exact_nt<ET> & a)
{
  ET e;
  internal::read_float_or_quotient(is, e);
  if (is)
    a = e;
  return is;
}

template< class ET >
class Is_valid< Lazy_exact_nt<ET> >
  : public std::unary_function< Lazy_exact_nt<ET>, bool > {
  public :
    bool operator()( const Lazy_exact_nt<ET>& x ) const {
      return is_valid(x.approx());
    }
};

template < typename ET >
struct NT_converter < Lazy_exact_nt<ET>, ET >
{
  const ET& operator()(const Lazy_exact_nt<ET> &a) const
  { return a.exact(); }
};

// Forward declaration to break inclusion cycle
namespace internal {
template<class>struct Exact_field_selector;
template<class>struct Exact_ring_selector;
}
// Compiler can deduce ET from the first argument.
template < typename ET >
struct NT_converter < Lazy_exact_nt<ET>,
  typename First_if_different<
    typename internal::Exact_field_selector<ET>::Type,
    ET>::Type>
{
  typename internal::Exact_field_selector<ET>::Type
    operator()(const Lazy_exact_nt<ET> &a) const
  { return NT_converter<ET,typename internal::Exact_field_selector<ET>::Type>()(a.exact()); }
};

template < typename ET >
struct NT_converter < Lazy_exact_nt<ET>,
  typename First_if_different<
   typename First_if_different<
    typename internal::Exact_ring_selector<ET>::Type,
    ET>::Type,
   typename internal::Exact_field_selector<ET>::Type>::Type>
{
  typename internal::Exact_ring_selector<ET>::Type operator()(const Lazy_exact_nt<ET> &a) const
  { return NT_converter<ET,typename internal::Exact_ring_selector<ET>::Type>()(a.exact()); }
};

namespace internal {
// Returns true if the value is representable by a double and to_double()
// would return it.  False means "don't know" (the exact number type is not
// queried).
template < typename ET >
inline bool
fit_in_double(const Lazy_exact_nt<ET>& l, double& r)
{ return fit_in_double(l.approx(), r); }

} // namespace internal

template <class NT_,class ROOT_, class ACDE_TAG_, class FP_TAG>
void
print(std::ostream &os, const CGAL::Lazy_exact_nt< Sqrt_extension<NT_,ROOT_,ACDE_TAG_,FP_TAG> > &r)
{
  print(os,r.exact());
}

namespace INTERN_LAZY_EXACT_NT {
template< typename ET , typename Tag>
class Modular_traits_base{
public:
  typedef Lazy_exact_nt<ET> NT;
  typedef ::CGAL::Tag_false Is_modularizable;
  typedef ::CGAL::Null_functor Residue_type;
  typedef ::CGAL::Null_functor Modular_image;  
  typedef ::CGAL::Null_functor Modular_image_representative;    
};

template< typename ET >
class Modular_traits_base<ET, Tag_true>{
  typedef Modular_traits<ET> MT_ET;
public:
  typedef Lazy_exact_nt<ET> NT;
  typedef CGAL::Tag_true Is_modularizable;
  typedef typename MT_ET::Residue_type Residue_type;

  struct Modular_image{
    Residue_type operator()(const NT& a){
      typename MT_ET::Modular_image modular_image;
      return modular_image(a.exact());
    }
  };
  struct Modular_image_representative{
    NT operator()(const Residue_type& x){
      typename MT_ET::Modular_image_representative modular_image_representative;
      return NT(modular_image_representative(x));
    }
  };    
};
} // namespace INTERN_LAZY_EXACT_NT

template < typename ET > 
class Modular_traits<Lazy_exact_nt<ET> >
  :public INTERN_LAZY_EXACT_NT::Modular_traits_base
<ET,typename Modular_traits<ET>::Is_modularizable>{};


#undef CGAL_double
#undef CGAL_int
#undef CGAL_To_interval

} //namespace CGAL

namespace Eigen {
  template<class> struct NumTraits;
  template<typename ET> struct NumTraits<CGAL::Lazy_exact_nt<ET> >
  {
    typedef CGAL::Lazy_exact_nt<ET> Real;
    // typedef CGAL::Lazy_exact_nt<ET> NonInteger;
    typedef CGAL::Lazy_exact_nt<typename NumTraits<ET>::NonInteger> NonInteger;
    typedef CGAL::Lazy_exact_nt<ET> Nested;
    typedef CGAL::Lazy_exact_nt<ET> Literal;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }

    enum {
      IsInteger = NumTraits<ET>::IsInteger,
      IsSigned = NumTraits<ET>::IsSigned,
      IsComplex = NumTraits<ET>::IsComplex,
      RequireInitialization = 1,
      ReadCost = 8,
      AddCost = 30,
      MulCost = 30
    };
  };
}

#endif // CGAL_LAZY_EXACT_NT_H
