// Copyright (c) 1999-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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


#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/number_utils.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/Number_type_traits.h>

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/utils.h>

#include <CGAL/Interval_nt.h>
#include <CGAL/Handle.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/NT_converter.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/Lazy_exact_nt_fwd.h>

#include <CGAL/Profile_counter.h>

#include <CGAL/functional_base.h> // Unary_function, Binary_function
#include <boost/iterator/transform_iterator.hpp> // for Root_of functor
#include <boost/static_assert.hpp>

#include <boost/operators.hpp>

#include <CGAL/Root_of_traits.h>

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
 * Lazy_exact_rep which has pure virtual methods .approx() and .exact().
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

CGAL_BEGIN_NAMESPACE

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


template <class ET>
class Lazy_exact_nt;

template <typename ET>
inline
void
print_dag(const Lazy_exact_nt<ET>& l, std::ostream& os, int level=0)
{
  l.print_dag(os, level);
}

inline
void
print_dag(double d, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++){
    os << "    ";
  }
  os << d << std::endl;
}


void
msg(std::ostream& os, int level, char* s)
  {
    int i;
    for(i = 0; i < level; i++){
      os << "    ";
    }
    os << s << std::endl;
  }

inline
void
print_dag(const Null_vector& nv, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++){
    os << "    ";
  }
  os << "Null_vector" << std::endl;
}

inline
void
print_dag(const Origin& nv, std::ostream& os, int level)
{
  for(int i = 0; i < level; i++){
    os << "    ";
  }
  os << "Origin" << std::endl;
}

#endif

// Abstract base class for lazy numbers and lazy objects
template <typename AT_, typename ET, typename E2A>
struct Lazy_construct_rep : public Rep
{
  typedef AT_ AT;

  AT at;
  mutable ET *et;

  Lazy_construct_rep ()
      : at(), et(NULL) {}

  Lazy_construct_rep (const AT& a)
      : at(a), et(NULL)
  {}

  Lazy_construct_rep (const AT& a, const ET& e)
      : at(a), et(new ET(e))
  {}

private:
  Lazy_construct_rep (const Lazy_construct_rep&) { std::abort(); } // cannot be copied.
public:

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
    }
  }

  virtual void print_dag(std::ostream& os, int level) const {}
#endif

  bool is_lazy() const { return et == NULL; }
  virtual void update_exact() = 0;
  virtual int depth() const  { return 1; }
  virtual ~Lazy_construct_rep () { delete et; };
};

// Abstract base representation class for lazy numbers
template <typename ET>
struct Lazy_exact_rep : public Lazy_construct_rep<Interval_nt<false>,
                                                  ET, CGAL_To_interval(ET) >
{
  typedef Lazy_construct_rep<Interval_nt<false>, ET, CGAL_To_interval(ET) > Base;

  Lazy_exact_rep (const Interval_nt<false> & i)
      : Base(i) {}

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
  }
#endif

private:
  Lazy_exact_rep (const Lazy_exact_rep&) { std::abort(); } // cannot be copied.

};

// int constant
template <typename ET>
struct Lazy_exact_Int_Cst : public Lazy_exact_rep<ET>
{
  Lazy_exact_Int_Cst (int i)
      : Lazy_exact_rep<ET>(double(i)) {}

  void update_exact()  { this->et = new ET((int)this->approx().inf()); }

};

// double constant
template <typename ET>
struct Lazy_exact_Cst : public Lazy_exact_rep<ET>
{
  Lazy_exact_Cst (double d)
      : Lazy_exact_rep<ET>(d) {}

  void update_exact()  { this->et = new ET(this->approx().inf()); }
};

// Exact constant
template <typename ET>
struct Lazy_exact_Ex_Cst : public Lazy_exact_rep<ET>
{
  Lazy_exact_Ex_Cst (const ET & e)
      : Lazy_exact_rep<ET>(CGAL_NTS to_interval(e))
  {
    this->et = new ET(e);
  }

  void update_exact()  { CGAL_assertion(false); }
};

// Construction from a Lazy_exact_nt<ET1> (which keeps the lazyness).
template <typename ET, typename ET1>
class Lazy_lazy_exact_Cst : public Lazy_exact_rep<ET>
{
  Lazy_exact_nt<ET1> l;

public:

  Lazy_lazy_exact_Cst (const Lazy_exact_nt<ET1> &x)
      : Lazy_exact_rep<ET>(x.approx()), l(x) {}

  void update_exact()
  {
    this->et = new ET(l.exact());
    this->approx() = l.approx();
    prune_dag();
  }
  int depth() const { return l.depth() + 1; }
  void prune_dag() { l = Lazy_exact_nt<ET1>::zero(); }
};


// Unary  operations: abs, sqrt, square.
// Binary operations: +, -, *, /, min, max.

// Base unary operation
template <typename ET>
struct Lazy_exact_unary : public Lazy_exact_rep<ET>
{
  Lazy_exact_nt<ET> op1;

  Lazy_exact_unary (const Interval_nt<false> &i, const Lazy_exact_nt<ET> &a)
      : Lazy_exact_rep<ET>(i), op1(a) {}

  int depth() const { return op1.depth() + 1; }
  void prune_dag() { op1 = Lazy_exact_nt<ET>::zero(); }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      msg(os, level, "Unary number operator:");
      print_dag(op1, os, level+1);
    }
  }
#endif
};

// Base binary operation
template <typename ET, typename ET1 = ET, typename ET2 = ET>
struct Lazy_exact_binary : public Lazy_exact_rep<ET>
{
  Lazy_exact_nt<ET1> op1;
  Lazy_exact_nt<ET2> op2;

  Lazy_exact_binary (const Interval_nt<false> &i,
		     const Lazy_exact_nt<ET1> &a, const Lazy_exact_nt<ET2> &b)
      : Lazy_exact_rep<ET>(i), op1(a), op2(b) {}

  int depth() const { return (std::max)(op1.depth(), op2.depth()) + 1; }
  void prune_dag()
  {
    op1 = Lazy_exact_nt<ET1>::zero();
    op2 = Lazy_exact_nt<ET2>::zero();
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      msg(os, level, "Binary number operator:");
      print_dag(op1, os, level+1);
      print_dag(op2, os, level+1);
    }
  }
#endif
};

// Here we could use a template class for all operations (STL provides
// function objects plus, minus, multiplies, divides...).  But it would require
// a template template parameter, and GCC 2.95 seems to crash easily with them.

#ifndef CGAL_CFG_COMMA_BUG
// Macro for unary operations
#define CGAL_LAZY_UNARY_OP(OP, NAME)                                     \
template <typename ET>                                                   \
struct NAME : public Lazy_exact_unary<ET>                                \
{                                                                        \
  typedef typename Lazy_exact_unary<ET>::AT::Protector P;                \
  NAME (const Lazy_exact_nt<ET> &a)                                      \
      : Lazy_exact_unary<ET>((P(), OP(a.approx())), a) {}                \
                                                                         \
  void update_exact()                                                    \
  {                                                                      \
    this->et = new ET(OP(this->op1.exact()));                            \
    if (!this->approx().is_point())                                      \
      this->approx() = CGAL_NTS to_interval(*(this->et));                \
    this->prune_dag();                                                   \
   }                                                                     \
};
#else
// Macro for unary operations
#define CGAL_LAZY_UNARY_OP(OP, NAME)                                     \
template <typename ET>                                                   \
struct NAME : public Lazy_exact_unary<ET>                                \
{                                                                        \
  typedef typename Lazy_exact_unary<ET>::AT::Protector P;                \
  NAME (const Lazy_exact_nt<ET> &a)                                      \
      : Lazy_exact_unary<ET>(a.approx() /* dummy value */, a)            \
  { P p; this->approx() = OP(a.approx()); }                              \
                                                                         \
  void update_exact()                                                    \
  {                                                                      \
    this->et = new ET(OP(this->op1.exact()));                            \
    if (!this->approx().is_point())                                      \
      this->approx() = CGAL_NTS to_interval(*(this->et));                \
    this->prune_dag();                                                   \
  }                                                                      \
};
#endif

CGAL_LAZY_UNARY_OP(opposite,  Lazy_exact_Opp)
CGAL_LAZY_UNARY_OP(CGAL_NTS abs,    Lazy_exact_Abs)
CGAL_LAZY_UNARY_OP(CGAL_NTS square, Lazy_exact_Square)
CGAL_LAZY_UNARY_OP(CGAL_NTS sqrt,   Lazy_exact_Sqrt)

#ifndef CGAL_CFG_COMMA_BUG
// A macro for +, -, * and /
#define CGAL_LAZY_BINARY_OP(OP, NAME)                                    \
template <typename ET, typename ET1 = ET, typename ET2 = ET>             \
struct NAME : public Lazy_exact_binary<ET, ET1, ET2>                     \
{                                                                        \
  typedef typename Lazy_exact_binary<ET,ET1,ET2>::AT::Protector P;	 \
  NAME (const Lazy_exact_nt<ET1> &a, const Lazy_exact_nt<ET2> &b)        \
    : Lazy_exact_binary<ET, ET1, ET2>((P(), a.approx() OP b.approx()), a, b) {} \
                                                                         \
  void update_exact()                                                    \
  {                                                                      \
    this->et = new ET(this->op1.exact() OP this->op2.exact());           \
    if (!this->approx().is_point())                                      \
      this->approx() = CGAL_NTS to_interval(*(this->et));                \
    this->prune_dag();                                                   \
   }                                                                     \
};
#else
// A macro for +, -, * and /
#define CGAL_LAZY_BINARY_OP(OP, NAME)                                    \
template <typename ET, typename ET1 = ET, typename ET2 = ET>             \
struct NAME : public Lazy_exact_binary<ET, ET1, ET2>                     \
{                                                                        \
  typedef typename Lazy_exact_binary<ET, ET1, ET2>::AT::Protector P;     \
  NAME (const Lazy_exact_nt<ET1> &a, const Lazy_exact_nt<ET2> &b)        \
    : Lazy_exact_binary<ET, ET1, ET2>(a.approx() /* dummy value */, a, b)\
  {P p; this->approx() = a.approx() OP b.approx(); }                     \
                                                                         \
  void update_exact()                                                    \
  {                                                                      \
    this->et = new ET(this->op1.exact() OP this->op2.exact());           \
    if (!this->approx().is_point())                                      \
      this->approx() = CGAL_NTS to_interval(*(this->et));                   \
    this->prune_dag();                                                   \
   }                                                                     \
};
#endif

CGAL_LAZY_BINARY_OP(+, Lazy_exact_Add)
CGAL_LAZY_BINARY_OP(-, Lazy_exact_Sub)
CGAL_LAZY_BINARY_OP(*, Lazy_exact_Mul)
CGAL_LAZY_BINARY_OP(/, Lazy_exact_Div)

// Minimum
template <typename ET>
struct Lazy_exact_Min : public Lazy_exact_binary<ET>
{
  Lazy_exact_Min (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_binary<ET>((min)(a.approx(), b.approx()), a, b) {}

  void update_exact()
  {
    this->et = new ET((min)(this->op1.exact(), this->op2.exact()));
    if (!this->approx().is_point()) this->approx() = CGAL_NTS to_interval(*(this->et));
    this->prune_dag();
  }
};

// Maximum
template <typename ET>
struct Lazy_exact_Max : public Lazy_exact_binary<ET>
{
  Lazy_exact_Max (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_binary<ET>((max)(a.approx(), b.approx()), a, b) {}

  void update_exact()
  {
    this->et = new ET((max)(this->op1.exact(), this->op2.exact()));
    if (!this->approx().is_point()) this->approx() = CGAL_NTS to_interval(*(this->et));
    this->prune_dag();
  }
};

// The real number type, handle class
template <typename ET_>
class Lazy_exact_nt
  : public Handle
  , boost::ordered_euclidian_ring_operators2< Lazy_exact_nt<ET_>, int >
  , boost::ordered_euclidian_ring_operators2< Lazy_exact_nt<ET_>, double >
{
public:

  typedef ET_                    ET;
  typedef Interval_nt_advanced   AT;

private:
  typedef Lazy_exact_nt<ET> Self;
  typedef Lazy_construct_rep< Interval_nt<false>, ET, To_interval<ET> > Self_rep;

public :

  typedef typename Number_type_traits<ET>::Has_gcd      Has_gcd;
  typedef typename Number_type_traits<ET>::Has_division Has_division;
  typedef typename Number_type_traits<ET>::Has_sqrt     Has_sqrt;

  typedef typename Number_type_traits<ET>::Has_exact_sqrt Has_exact_sqrt;
  typedef typename Number_type_traits<ET>::Has_exact_division
                                                        Has_exact_division;
  typedef typename Number_type_traits<ET>::Has_exact_ring_operations
                                                     Has_exact_ring_operations;

  Lazy_exact_nt (Self_rep *r)
  { PTR = r; }

  Lazy_exact_nt ()
    : Handle(zero()) {}

  Lazy_exact_nt (const CGAL_int(ET) & i)
  { PTR = new Lazy_exact_Int_Cst<ET>(i); }

  Lazy_exact_nt (unsigned i)
  { PTR = new Lazy_exact_Cst<ET>(i); }

  Lazy_exact_nt (const CGAL_double(ET) & d)
  { PTR = new Lazy_exact_Cst<ET>(d); }

  Lazy_exact_nt (const ET & e)
  { PTR = new Lazy_exact_Ex_Cst<ET>(e); }

  template <class ET1>
  Lazy_exact_nt (const Lazy_exact_nt<ET1> &x)
  { PTR = new Lazy_lazy_exact_Cst<ET, ET1>(x); }

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
    ET res = exact();
    res %= b.exact();
    return *this = Lazy_exact_nt<ET>(res);
  }

  Self & operator%=(int b)
  {
    CGAL_precondition(b != 0);
    ET res = exact();
    res %= b;
    return *this = Lazy_exact_nt<ET>(res);
  }

  Interval_nt<true> interval() const
  {
    const Interval_nt<false>& i = approx();
    return Interval_nt<true>(i.inf(), i.sup());
  }

  const Interval_nt<false>& approx() const
  { return ptr()->approx(); }

  Interval_nt_advanced approx_adv() const
  { return ptr()->approx(); }

  const ET & exact() const
  { return ptr()->exact(); }

  int depth() const
  { return ptr()->depth(); }

  void
  print_dag(std::ostream& os, int level) const
  {
    ptr()->print_dag(os, level);
  }

  static const double & get_relative_precision_of_to_double()
  {
      return relative_precision_of_to_double;
  }

  static void set_relative_precision_of_to_double(const double & d)
  {
      CGAL_assertion(d > 0 && d < 1);
      relative_precision_of_to_double = d;
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

  // We have a static variable for optimizing zero and default constructor.
  static const Self & zero()
  {
    static const Self z = new Lazy_exact_Int_Cst<ET>(0);
    return z;
  }

private:
  Self_rep * ptr() const { return (Self_rep*) PTR; }

  static double relative_precision_of_to_double;
};


template <typename ET>
double Lazy_exact_nt<ET>::relative_precision_of_to_double = 0.00001;


template <typename ET1, typename ET2>
bool
operator<(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  if (a.identical(b))
    return false;
  Uncertain<bool> res = a.approx() < b.approx();
  if (is_singleton(res))
    return res;
  CGAL_PROFILER(std::string("failures of : ") + std::string(CGAL_PRETTY_FUNCTION));
  return a.exact() < b.exact();
}

template <typename ET1, typename ET2>
bool
operator==(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  if (a.identical(b))
    return true;
  Uncertain<bool> res = a.approx() == b.approx();
  if (is_singleton(res))
    return res;
  CGAL_PROFILER(std::string("failures of : ") + std::string(CGAL_PRETTY_FUNCTION));
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
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  Uncertain<bool> res = a.approx() < b;
  if (is_singleton(res))
    return res;
  CGAL_PROFILER(std::string("failures of : ") + std::string(CGAL_PRETTY_FUNCTION));
  return a.exact() < b;
}

template <typename ET>
bool
operator>(const Lazy_exact_nt<ET>& a, int b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  Uncertain<bool> res = b < a.approx();
  if (is_singleton(res))
    return res;
  CGAL_PROFILER(std::string("failures of : ") + std::string(CGAL_PRETTY_FUNCTION));
  return b < a.exact();
}

template <typename ET>
bool
operator==(const Lazy_exact_nt<ET>& a, int b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  Uncertain<bool> res = b == a.approx();
  if (is_singleton(res))
    return res;
  CGAL_PROFILER(std::string("failures of : ") + std::string(CGAL_PRETTY_FUNCTION));
  return b == a.exact();
}


template <typename ET1, typename ET2>
Lazy_exact_nt< typename Coercion_traits<ET1, ET2>::Coercion_type >
operator+(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  return new Lazy_exact_Add<typename Coercion_traits<ET1, ET2>::Coercion_type,
                            ET1, ET2>(a, b);
}

template <typename ET1, typename ET2>
Lazy_exact_nt< typename Coercion_traits<ET1, ET2>::Coercion_type >
operator-(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  return new Lazy_exact_Sub<typename Coercion_traits<ET1, ET2>::Coercion_type,
                            ET1, ET2>(a, b);
}

template <typename ET1, typename ET2>
Lazy_exact_nt< typename Coercion_traits<ET1, ET2>::Coercion_type >
operator*(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  return new Lazy_exact_Mul<typename Coercion_traits<ET1, ET2>::Coercion_type,
                            ET1, ET2>(a, b);
}

template <typename ET1, typename ET2>
Lazy_exact_nt< typename Coercion_traits<ET1, ET2>::Coercion_type >
operator/(const Lazy_exact_nt<ET1>& a, const Lazy_exact_nt<ET2>& b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  CGAL_precondition(b != 0);
  return new Lazy_exact_Div<typename Coercion_traits<ET1, ET2>::Coercion_type,
                            ET1, ET2>(a, b);
}

//
// Algebraic structure traits
//


namespace INTERN_LAZY_EXACT_NT {
  
template< class NT, class Functor >
struct Simplify_selector {
  struct Simplify : public Unary_function<NT&, void> {
    void operator()( NT& x ) const {
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
  struct Unit_part : public Unary_function<NT, NT > {
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
  struct Is_zero : public Unary_function<NT, bool > {
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
  struct Is_one : public Unary_function<NT, bool > {
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
  struct Square : public Unary_function<NT, NT > {  
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
  struct Integral_division : public Binary_function<NT, NT, NT > {
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
  struct Is_square : public Binary_function<NT, NT&, bool > {
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
    struct Sqrt : public Unary_function<NT, NT > { 
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
  struct Kth_root : public Binary_function<int, NT, NT > {
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
      typedef Arity_tag< 3 >         Arity;
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
  struct Gcd : public Binary_function<NT, NT, NT > {
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
  struct Div : public Binary_function<NT, NT, NT > {
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
struct Mod_selector {
  struct Mod : public Binary_function<NT, NT, NT > {
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
    typedef Arity_tag< 4 >         Arity;
    
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
      BOOST_STATIC_ASSERT((::boost::is_same<
        typename Coercion_traits< NT1, NT2 >::Coercion_type, NT
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
  
  
  
/*  template< class ET, class Algebraic_structure_tag >
  class Lazy_exact_algebraic_structure_traits_base;
  
  template< class ET, class Algebraic_structure_tag >
  struct Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>, 
                                                    Algebraic_structure_tag >
    : public Algebraic_structure_traits_base< Lazy_exact_nt<ET>, 
                                              Algebraic_structure_tag > {
        typedef Lazy_exact_nt<ET> Algebraic_structure;
        typedef typename Algebraic_structure_traits<ET>::Is_exact Is_exact; 
        
      class Square 
        : public Unary_function< Algebraic_structure, Algebraic_structure > {
        public:
          Algebraic_structure operator()( const Algebraic_structure& x ) const {
              CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
              return new Lazy_exact_Square<ET>(x);
          }
      };    
  };

  template< class ET >
  class Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>,
                                 Integral_domain_without_division_tag >
    : public Algebraic_structure_traits_base< Lazy_exact_nt<ET>, 
                                 Integral_domain_without_division_tag > {
    public:
      typedef typename Algebraic_structure_traits<ET>::Is_exact Is_exact; 
      typedef Lazy_exact_nt<ET> Algebraic_structure;
      class Unit_part
        : public Unary_function< Algebraic_structure, Algebraic_structure > {
        public:
          Algebraic_structure operator()( const Algebraic_structure& x ) const {
            return Lazy_exact_nt<ET>( CGAL_NTS unit_part( x.exact() ) );
          }
      };
  
      class Is_zero
        : public Unary_function< Algebraic_structure, bool > {
        public:
          bool operator()( const Algebraic_structure& x ) const {
            return CGAL_NTS is_zero( x.exact() );
          }
      };
  
      class Is_one
        : public Unary_function< Algebraic_structure, bool > {
        public:
          bool operator()( const Algebraic_structure& x ) const {
            return CGAL_NTS is_one( x.exact() );
          }
      };
  
      class Square 
        : public Unary_function< Algebraic_structure, Algebraic_structure > {
        public:
          Algebraic_structure operator()( const Algebraic_structure& x ) const {
              CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
              return new Lazy_exact_Square<ET>(x);
          }
      };    
  };
  
  template< class ET >
  class Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>,
                                                    Integral_domain_tag >
    : public Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>, 
                                  Integral_domain_without_division_tag > {
    public:
      typedef Lazy_exact_nt<ET>          Algebraic_structure;
      typedef Integral_domain_tag  Algebraic_structure_tag;
      
      class Integral_division
        : public Binary_function< Algebraic_structure, Algebraic_structure,
                                  Algebraic_structure > {
        public:
          Algebraic_structure operator()( const Algebraic_structure& x,
                                          const Algebraic_structure& y ) const {
            return Lazy_exact_nt<ET>( CGAL_NTS integral_division( x.exact(), 
                                                                  y.exact() ) );
          }
          
          CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
      };
  };

  template< class ET >
  class Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>,
                                                    Field_tag >
    : public Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>, 
                                                  Integral_domain_tag > {
    public:
      typedef Lazy_exact_nt<ET>  Algebraic_structure;
      typedef Field_tag    Algebraic_structure_tag;
  };

  template< class ET >
  class Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>,
                                                    Field_with_sqrt_tag >
    : public Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>, 
                                                         Field_tag > {
    public:
      typedef Lazy_exact_nt<ET>          Algebraic_structure;
      typedef Field_with_sqrt_tag  Algebraic_structure_tag;
                                                                                   
      class Sqrt 
        : public Unary_function< Algebraic_structure, Algebraic_structure > {
        public:
          Algebraic_structure operator()( const Algebraic_structure& x ) const {
            CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
            CGAL_precondition(x >= 0);
            return new Lazy_exact_Sqrt<ET>(x);
          }
      };
      
  };

  template< class ET >
  class Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>,
                                                Field_with_kth_root_tag >
    : public Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>, 
                                                  Field_with_sqrt_tag > {

    public:
      typedef Lazy_exact_nt<ET>          Algebraic_structure;
      typedef Field_with_kth_root_tag  Algebraic_structure_tag;
                                                                                         
      class Kth_root
        : public Binary_function< int, Algebraic_structure, 
                                  Algebraic_structure > {
        public:
          Algebraic_structure operator()( int k, 
                                          const Algebraic_structure& x ) const {
            return Lazy_exact_nt<ET>( CGAL_NTS kth_root( k, x.exact() ) );
          }
      };
      
  };

  template< class ET >
  class Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>,
                                                Field_with_root_of_tag >
    : public Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>, 
                                              Field_with_kth_root_tag > {
    private:
      struct Cast{                                      
        typedef ET result_type;                               
        ET operator()(const Lazy_exact_nt<ET>& lazy_exact) const { 
          return lazy_exact.exact();
        }
      }; 

    public:
      typedef Lazy_exact_nt<ET>          Algebraic_structure;
      typedef Field_with_root_of_tag  Algebraic_structure_tag;
                                                                                         
      class Root_of {
        public:
          typedef Algebraic_structure result_type;
          typedef Arity_tag< 3 >         Arity;

          template< class Input_iterator >
          Algebraic_structure operator()( int k, Input_iterator begin,
                                                 Input_iterator end ) const {
            Cast cast;
            return Lazy_exact_nt<ET>( typename Algebraic_structure_traits<ET>::
                                                  Root_of()( k, 
                              ::boost::make_transform_iterator( begin, cast ), 
                              ::boost::make_transform_iterator( end, cast ) ) );
          }
      };      
  };

  template< class ET >
  class Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>,
                                        Unique_factorization_domain_tag >
    : public Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>, 
                                              Integral_domain_tag > {
    public:
      typedef Lazy_exact_nt<ET>          Algebraic_structure;
      typedef Unique_factorization_domain_tag  Algebraic_structure_tag;

      // gcd kills filtering.
      class Gcd 
        : public Binary_function< Algebraic_structure, Algebraic_structure,
                                  Algebraic_structure > {
        public:
          Algebraic_structure operator()( const Algebraic_structure& x,
                                         const Algebraic_structure& y ) const {
              CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
              return Lazy_exact_nt<ET>( CGAL_NTS gcd( x.exact(), y.exact() ) );
          }
          
          CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
      };                                                                                         
  };

  template< class ET >
  class Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>,
                                                    Euclidean_ring_tag >
    : public Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>, 
                                      Unique_factorization_domain_tag > {
    public:
      typedef Lazy_exact_nt<ET>          Algebraic_structure;
      typedef Euclidean_ring_tag  Algebraic_structure_tag;

      class Div
        : public Binary_function< Algebraic_structure, Algebraic_structure,
                                  Algebraic_structure > {
        public:
          Algebraic_structure operator()( const Algebraic_structure& x,
                                         const Algebraic_structure& y ) const {
            return Lazy_exact_nt<ET>( CGAL_NTS div( x.exact(), y.exact() ) );
          }
          
          CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
      };                                                                                         

      class Mod
        : public Binary_function< Algebraic_structure, Algebraic_structure,
                                  Algebraic_structure > {
        public:
          Algebraic_structure operator()( const Algebraic_structure& x,
                                         const Algebraic_structure& y ) const {
            return Lazy_exact_nt<ET>( CGAL_NTS mod( x.exact(), y.exact() ) );
          }
          
          CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Algebraic_structure )
      };
      
      class Div_mod {
        public:
          typedef void                  result_type;
          typedef Algebraic_structure   first_argument_type;
          typedef Algebraic_structure   second_argument_type;
          typedef Algebraic_structure&  third_arugment_type;
          typedef Algebraic_structure&  fourth_argument_type;
          typedef Arity_tag< 4 >         Arity;
          
          void operator()( const Algebraic_structure& x, 
                           const Algebraic_structure& y,
                           Algebraic_structure& q, 
                           Algebraic_structure& r ) const {
            ET q_et;
            ET r_et;
            CGAL_NTS div_mod( x.exact(), y.exact(), q_et, r_et );
            q = Lazy_exact_nt<ET>( q_et );
            r = Lazy_exact_nt<ET>( r_et );
          }
          
          template< class NT1, class NT2 >
          void operator()( const NT1& x, const NT2& y,
                           Algebraic_structure& q, 
                           Algebraic_structure& r ) const {
            BOOST_STATIC_ASSERT((::boost::is_same<
                 typename Coercion_traits< NT1, NT2 >::Coercion_type, Algebraic_structure
                                                 >::value));
            
            typename Coercion_traits< NT1, NT2 >::Cast cast;
            operator()( cast(x), cast(y), q, r );                      
          }
      };
  };*/  
} // INTERN_LAZY_EXACT_NT

     

/*template< class ET > class Algebraic_structure_traits< Lazy_exact_nt<ET> >
  : public INTERN_LAZY_EXACT_NT::
                Lazy_exact_algebraic_structure_traits_base< Lazy_exact_nt<ET>,
         typename Algebraic_structure_traits<ET>::Algebraic_structure_tag > {
};*/


template <class ET>
class Algebraic_structure_traits< Lazy_exact_nt<ET> >
    :public Algebraic_structure_traits_base
      < Lazy_exact_nt<ET>, 
       typename Algebraic_structure_traits<ET>::Algebraic_structure_tag >
{
private:
    typedef Algebraic_structure_traits<ET> AST_ET;
    typedef typename AST_ET::Algebraic_structure_tag ET_as_tag;
public:
    typedef typename AST_ET::Is_exact Is_exact;

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
};



//
// Real embeddalbe traits
//

template < typename ET > class Real_embeddable_traits< Lazy_exact_nt<ET> > 
  : public Real_embeddable_traits_base< Lazy_exact_nt<ET> > {
  
  // Every type ET of Lazy_exact_nt<ET> has to be real embeddable.
  BOOST_STATIC_ASSERT((::boost::is_same< typename Real_embeddable_traits< ET >
                                ::Is_real_embeddable, Tag_true >::value));

  public:
    typedef Lazy_exact_nt<ET> Real_embeddable;  
           
    class Abs 
      : public Unary_function< Real_embeddable, Real_embeddable > {
      public:
        Real_embeddable operator()( const Real_embeddable& a ) const {
            CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
            return new Lazy_exact_Abs<ET>(a);
        }
    };
    
    class Sign 
      : public Unary_function< Real_embeddable, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Real_embeddable& a ) const {
            CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
            Uncertain< ::CGAL::Sign> res = CGAL_NTS sign(a.approx());
            if (is_singleton(res))
                return res;
            CGAL_PROFILER(std::string("failures of : ") + std::string(CGAL_PRETTY_FUNCTION));
            return CGAL_NTS sign(a.exact());
        }
    };
    
    class Compare 
      : public Binary_function< Real_embeddable, Real_embeddable,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Real_embeddable& a, 
                                            const Real_embeddable& b ) const {
            CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
            if (a.identical(b))
                return EQUAL;
            Uncertain<Comparison_result> res = CGAL_NTS compare(a.approx(), b.approx());
            if (is_singleton(res))
                return res;
            CGAL_PROFILER(std::string("failures of : ") + std::string(CGAL_PRETTY_FUNCTION));
            return CGAL_NTS compare(a.exact(), b.exact());
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Real_embeddable,
                                                      Comparison_result );
        
    };
    
    class To_double 
      : public Unary_function< Real_embeddable, double > {
      public:
        double operator()( const Real_embeddable& a ) const {
            CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
            
            const Interval_nt<false>& app = a.approx();
            if (app.sup() == app.inf())
                return app.sup();
            
            // If it's precise enough, then OK.
            if ((app.sup() - app.inf())
                    < Lazy_exact_nt<ET>::get_relative_precision_of_to_double()
                    * (std::max)(std::fabs(app.inf()), std::fabs(app.sup())) )
                return CGAL_NTS to_double(app);
            
            CGAL_PROFILER(std::string("failures of : ") + std::string(CGAL_PRETTY_FUNCTION));
            
            // Otherwise we trigger exact computation first,
            // which will refine the approximation.
            a.exact();
            return CGAL_NTS to_double(a.approx());
        }
    };
    
    class To_interval 
      : public Unary_function< Real_embeddable, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Real_embeddable& a ) const {
            CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
            return a.approx().pair();
        }
    };
    
    class Is_finite
      : public Unary_function< Real_embeddable, bool > {
      public:
        bool operator()( const Real_embeddable& x ) const {
          return CGAL_NTS is_finite(x.approx()) || CGAL_NTS is_finite(x.exact());        
        }
    };
    
};

template <class ET1, class ET2, class F>
class Lazy_exact_nt_coercion_traits_base{
    typedef Tag_false Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
    //typedef Null_type    Coercion_type
    typedef Null_functor Cast;
};

template <class ET1, class ET2>
class Lazy_exact_nt_coercion_traits_base
< Lazy_exact_nt<ET1>, Lazy_exact_nt<ET2>, Tag_true >{
    typedef Coercion_traits<ET1,ET2> CT;
    typedef Lazy_exact_nt<ET1> A;
    typedef Lazy_exact_nt<ET2> B;
public:
    typedef Lazy_exact_nt<typename CT::Coercion_type> Coercion_type;
    typedef typename CT::Are_implicit_interoperable Are_explicit_interoperable; 
    typedef typename CT::Are_implicit_interoperable Are_implicit_interoperable; 
    
    class Cast{
    private:
        template <class T>
        inline Coercion_type cast(const T& x) const{ return Coercion_type(x); }
        inline Coercion_type cast(const Coercion_type& x) const{ return x; }
    public:
        typedef Coercion_type result_type;       
        Coercion_type operator()(const A& x) const { return cast(x);}
        Coercion_type operator()(const B& x) const { return cast(x);}         
    };
};


CGAL_DEFINE_COERCION_TRAITS_FOR_SELF_TEM(Lazy_exact_nt<ET>, class ET);

template<class ET1, class ET2 >
class Coercion_traits< Lazy_exact_nt<ET1>, Lazy_exact_nt<ET2> >
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
        ::type  Coercion_type;                                          \
        typedef typename boost::mpl::if_c <interoperable, Null_functor, \
    INTERN_CT::Cast_from_to<NTX,NT> >::type Cast;                       \
    };                                                                  \
                                                                        \
    template<class ET>                                                  \
    struct Coercion_traits< Lazy_exact_nt<ET>, NTX >                    \
        :public Coercion_traits<NTX, Lazy_exact_nt<ET> >{};             \
    

CGAL_COERCION_TRAITS_LAZY_EXACT(int);
CGAL_COERCION_TRAITS_LAZY_EXACT(short);
CGAL_COERCION_TRAITS_LAZY_EXACT(double);
CGAL_COERCION_TRAITS_LAZY_EXACT(float);
#undef CGAL_COERCION_TRAITS_LAZY_EXACT

template <typename ET>
inline
Lazy_exact_nt<ET>
min BOOST_PREVENT_MACRO_SUBSTITUTION (const Lazy_exact_nt<ET> & a,
				      const Lazy_exact_nt<ET> & b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  return new Lazy_exact_Min<ET>(a, b);
}

template <typename ET>
inline
Lazy_exact_nt<ET>
max BOOST_PREVENT_MACRO_SUBSTITUTION (const Lazy_exact_nt<ET> & a,
				      const Lazy_exact_nt<ET> & b)
{
  CGAL_PROFILER(std::string("calls to    : ") + std::string(CGAL_PRETTY_FUNCTION));
  return new Lazy_exact_Max<ET>(a, b);
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
  is >> e;
  if (is)
    a = e;
  return is;
}

template< class ET >
class Is_valid< Lazy_exact_nt<ET> > 
  : public Unary_function< Lazy_exact_nt<ET>, bool > {
  public :
    bool operator()( const Lazy_exact_nt<ET>& x ) const {
      return is_valid(x.approx()) || is_valid(x.exact());
    }  
};

template <typename ET>
inline
io_Operator
io_tag (const Lazy_exact_nt<ET>&)
{ return io_Operator(); }

template < typename ET >
struct NT_converter < Lazy_exact_nt<ET>, ET >
{
  const ET& operator()(const Lazy_exact_nt<ET> &a) const
  { return a.exact(); }
};

// Returns true if the value is representable by a double and to_double()
// would return it.  False means "don't know".
template < typename ET >
inline bool
fit_in_double(const Lazy_exact_nt<ET>& l, double& r)
{ return fit_in_double(l.approx(), r); }


// We create a type of new node in Lazy_exact_nt's DAG
// for the make_root_of_2() operation.

#if 0 // To be finished
template <typename ET >
struct Lazy_exact_ro2
  : public Lazy_exact_rep< typename Root_of_traits<ET>::RootOf_2 >
{
    typedef typename Root_of_traits<ET>::RootOf_2   RO2;
    typedef Lazy_exact_rep<RO2>                     Base;
    typedef typename Base::AT::Protector            P;


    mutable Lazy_exact_nt<ET> op1, op2, op3;
    bool smaller;

    Lazy_exact_ro2 (const Lazy_exact_nt<ET> &a,
                    const Lazy_exact_nt<ET> &b,
                    const Lazy_exact_nt<ET> &c, bool s)
#ifndef CGAL_CFG_COMMA_BUG
      : Base((P(), make_root_of_2(a.approx(), b.approx(), c.approx(), s))),
        op1(a), op2(b), op3(c), smaller(s) {}
#else
      : Base(a.approx() /* dummy value */, a),
        op1(a), op2(b), op3(c), smaller(s)
  {
    P p;
    this->approx() = make_root_of_2(a.approx(), b.approx(),
                                    c.approx(), s);
  }
#endif

    void update_exact()
    {
        this->et = new RO2(make_root_of_2(op1.exact(), op2.exact(),
                                          op3.exact(), smaller));

        if (!this->approx().is_point())
            this->at = to_interval(*(this->et));
        this->prune_dag();

    }

    void prune_dag() const
    {
        op1 = op2 = op3 = Lazy_exact_nt<ET>::zero();
    }
};

template < typename ET >
inline
Lazy_exact_nt< typename Root_of_traits<ET>::RootOf_2 >
make_root_of_2( const Lazy_exact_nt<ET> &a,
                const Lazy_exact_nt<ET> &b,
                const Lazy_exact_nt<ET> &c, bool d)
{
    return new Lazy_exact_ro2<ET>(a, b, c, d);
}

template <typename NT >
struct Root_of_traits< Lazy_exact_nt < NT > >
{
private:
    typedef Root_of_traits<NT> T;
public:
    typedef Lazy_exact_nt< typename T::RootOf_1 > RootOf_1;
    typedef Lazy_exact_nt< typename T::RootOf_2 > RootOf_2;
};

#endif // 0

#undef CGAL_double
#undef CGAL_int
#undef CGAL_To_interval

CGAL_END_NAMESPACE

#endif // CGAL_LAZY_EXACT_NT_H
