// ============================================================================
//
// Copyright (c) 1999,2000 The CGAL Consortium
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
// file          : include/CGAL/Lazy_exact_nt.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_LAZY_EXACT_NT_H
#define CGAL_LAZY_EXACT_NT_H

#include <CGAL/basic.h>
#include <CGAL/number_utils.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Handle_for.h>

/*
 * This file contains the definition of the number type Lazy_exact_nt<ET>,
 * where ET is an exact number type (must provide the exact operations needed).
 *
 * Lazy_exact_nt<ET> provides a DAG-based lazy evaluation, like LEDA's real,
 * Core's Expr, LEA's lazy rationals...
 *
 * The values are first approximated using Interval_nt.
 * The exactness is provided when needed by ET.
 *
 * Lazy_exact_nt<ET> is just a handle to the abstract base class
 * Lazy_exact_nt_rep which has pure virtual methods .approx() and .exact().
 * From this class derives one class per operation, with one constructor.
 *
 * The DAG is managed by :
 * - Handle_for<RefCounted, Allocator> and RefCounted.
 * - virtual functions to denote the various operators (instead of an enum).
 *
 * Other packages with vaguely similar design : APU, MetaCGAL, LOOK.
 */

/*
 * TODO (vaguely by decreasing priority):
 * - A choice must be made between Interval_nt/Interval_nt_advanced (or both)
 *   for access members.
 * - There's probably a mean to factorize something in the base class for
 *   .exact().
 ------ [ passed this, it'll be usable ] -----
 * - Don't use the CGAL template for sign, compare.
 * - Do we want to have an interval rafinement functionnality ?
 * - The next step will be to replace Cartesian<Lazy_exact_nt<ET> > by a kernel
 *   similar to LOOK in functionality : Lazy_exact_Cartesian<FT,ET>.
 * - Predicates could use the filtered advanced version.
 * - Geometric constructions could use the interval_advanced.
 ------ [ style issues ] -----
 * - Make all binary/unary operations via a template, having the template
 *   operation as argument (needs a function objet for that...).
 * - Separate the handle and the representation in 2 files.
 * - Add an Allocator template argument ?
 * - Add a CT template parameter like Filtered_exact_nt<> ?
 * - Add a string constant to provide an expression string (a la MetaCGAL) ?
 *   // virtual ostream operator<<() const = 0; // or string, like Core ?
 */

/*
 * Choices that have to be made about .approx():
 * - virtual or not ?
 * - lazy evaluation or not ?
 *
 *                        lazy                  no-lazy
 * virtual                                    
 * no-virtual            impossible           
 *
 * lazy        => virtual
 * lazy        => less rounding mode changes
 * virtual     => less space for constants (for double we loose only 8bytes)
 *                more flexibility
 *             => similar implementation with exact()
 *             => allows rafinement schemes
 *
 * non virtual => more straightforward to code ?
 * lazy        => need to use a boolean "status" or an invalid interval.
 *
 * I think the question is what costs more between:
 * - a virtual function call, and
 * - changing the rounding mode.
 *
 * For the moment, I choose to implement it lazily (=> virtual) ???
 */

CGAL_BEGIN_NAMESPACE

template <typename ET> class Lazy_exact_nt;

// Abstract base representation class
template <typename ET>
struct Lazy_exact_nt_rep : public Ref_counted
{
  mutable Interval_nt_advanced in;

  Lazy_exact_nt_rep ()
	  : in(1,0) {}
  Lazy_exact_nt_rep (const double d)
	  : in(d) {}
  // virtual void update_interval() const = 0;
  // update_interval() is supposed to be called with rounding -> +inf.
  virtual void update_approx() const = 0;
  virtual ET exact() const = 0;
  Interval_nt_advanced interval() const
  {
      if (!is_valid(in))
      {
	  FPU_CW_t backup = FPU_get_and_set_cw(CGAL_FE_UPWARD);
	  update_approx();
	  FPU_set_cw(backup);
      }
      return in;
  }
  Interval_nt_advanced approx() const
  {
      if (!is_valid(in))
	  update_approx();
      return in;
  }
  virtual ~Lazy_exact_nt_rep () {}; // ok ?
};

// double constant
template <typename ET>
struct Lazy_exact_nt_Cst : public Lazy_exact_nt_rep<ET>
{
  Lazy_exact_nt_Cst (const double a)
      : Lazy_exact_nt_rep<ET>(a) {}

  void update_approx() const
  {
	  CGAL_assertion(false);
  }
  ET exact() const
  {
	  return ET(in.inf());
  }
  // ~Lazy_exact_nt_Cst() {}
};

// Unary operations: abs, sqrt, square.

// Base unary operation
template <typename ET>
struct Lazy_exact_nt_unary : public Lazy_exact_nt_rep<ET>
{
  mutable ET *et;
  const Lazy_exact_nt<ET> op1;

  Lazy_exact_nt_unary (const Lazy_exact_nt<ET> &a)
      : Lazy_exact_nt_rep<ET>(), et(NULL), op1(a) {}

  // Interval_nt_advanced interval() const {};
};

// abs
template <typename ET>
struct Lazy_exact_nt_Abs : public Lazy_exact_nt_unary<ET>
{
  Lazy_exact_nt_Abs (const Lazy_exact_nt<ET> &a)
      : Lazy_exact_nt_unary<ET>(a) {}

  void update_approx() const
  {
      in = CGAL_NTS abs(op1.approx());
  }
  ET exact() const
  {
      if (!et)
	  et = new ET(CGAL_NTS abs(op1.exact()));
      return *et;
  }
};

// sqrt
template <typename ET>
struct Lazy_exact_nt_Sqrt : public Lazy_exact_nt_unary<ET>
{
  Lazy_exact_nt_Sqrt (const Lazy_exact_nt<ET> &a)
      : Lazy_exact_nt_unary<ET>(a) {}

  void update_approx() const
  {
      in = sqrt(op1.approx());
  }
  ET exact() const
  {
      if (!et)
	  et = new ET(CGAL::sqrt(op1.exact()));
      return *et;
  }
};

// square
template <typename ET>
struct Lazy_exact_nt_Square : public Lazy_exact_nt_unary<ET>
{
  Lazy_exact_nt_Square (const Lazy_exact_nt<ET> &a)
      : Lazy_exact_nt_unary<ET>(a) {}

  void update_approx() const
  {
      in = CGAL_NTS square(op1.approx());
  }
  ET exact() const
  {
      if (!et)
	  et = new ET(CGAL_NTS square(op1.exact()));
      return *et;
  }
};


// Binary operations: +, -, *, /, min, max.

// Base binary operation (should it derive from unary instead ?)
template <typename ET>
struct Lazy_exact_nt_binary : public Lazy_exact_nt_rep<ET>
{
  mutable ET *et;
  const Lazy_exact_nt<ET> op1, op2;
  Lazy_exact_nt_binary (const Lazy_exact_nt<ET> &a,
	                    const Lazy_exact_nt<ET> &b)
      : Lazy_exact_nt_rep<ET>(), et(NULL), op1(a), op2(b) {}
  // ~Lazy_exact_nt_binary() { }
};

#if 0
// Template binary operator (might be merged with the above ?)
// Note : G++ 2.95 produces an ICE on this  :(  Will try again later...
template <typename ET, template <typename T> class Op>
struct Lazy_exact_nt_binary_op : public Lazy_exact_nt_binary<ET>
{
  Lazy_exact_nt_binary_op (const Lazy_exact_nt<ET> &a,
		               const Lazy_exact_nt<ET> &b)
    : Lazy_exact_nt_binary<ET>(a,b) {}

  void update_approx() const
  {
      in = Op<Interval_nt_advanced>()(op1.approx(), op2.approx());
  }
  ET exact() const
  {
      if (!et)
	  et = new ET(Op<ET>()(op1.exact(), op2.exact()));
      return *et;
  }
};

#elsif 0
// Second solution : the macro :)

#define CGAL_LAZY_BINARY_OP(OP)                                       \
template <typename ET>                                                \
struct Lazy_exact_nt_##OP : public Lazy_exact_nt_binary<ET>   \
{                                                                     \
  Lazy_exact_nt_##OP (const Lazy_exact_nt<ET> &a,                 \
		          const Lazy_exact_nt<ET> &b)                 \
    : Lazy_exact_nt_binary<ET>(a,b) {}                            \
                                                                      \
  void update_approx() const                                          \
  {                                                                   \
      in = OP<Interval_nt_advanced>()(op1.approx(), op2.approx());    \
  }                                                                   \
  ET exact() const                                                    \
  {                                                                   \
      if (!et)                                                        \
	  et = new ET(OP<ET>()(op1.exact(), op2.exact()));            \
      return *et;                                                     \
  }                                                                   \
};

CGAL_LAZY_BINARY_OP(Min)

#endif

// Addition
template <typename ET>
struct Lazy_exact_nt_Add : public Lazy_exact_nt_binary<ET>
{
  Lazy_exact_nt_Add (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_nt_binary<ET>(a,b) {}

  void update_approx() const
  {
      in = op1.approx() + op2.approx();
  }
  ET exact() const
  {
      if (!et)
	  et = new ET(op1.exact() + op2.exact());
      return *et;
  }
};

// Subtraction
template <typename ET>
struct Lazy_exact_nt_Sub : public Lazy_exact_nt_binary<ET>
{
  Lazy_exact_nt_Sub (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_nt_binary<ET>(a,b) {}

  void update_approx() const
  {
      in = op1.approx() - op2.approx();
  }
  ET exact() const
  {
      if (!et)
	  et = new ET(op1.exact() - op2.exact());
      return *et;
  }
};

// Multiplication
template <typename ET>
struct Lazy_exact_nt_Mul : public Lazy_exact_nt_binary<ET>
{
  Lazy_exact_nt_Mul (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_nt_binary<ET>(a,b) {}

  void update_approx() const
  {
      in = op1.approx() * op2.approx();
  }
  ET exact() const
  {
      if (!et)
	  et = new ET(op1.exact() * op2.exact());
      return *et;
  }
};

// Division
template <typename ET>
struct Lazy_exact_nt_Div : public Lazy_exact_nt_binary<ET>
{
  Lazy_exact_nt_Div (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_nt_binary<ET>(a,b) {}

  void update_approx() const
  {
      in = op1.approx() / op2.approx();
  }
  ET exact() const
  {
      if (!et)
	  et = new ET(op1.exact() / op2.exact());
      return *et;
  }
};

// Minimum
template <typename ET>
struct Lazy_exact_nt_Min : public Lazy_exact_nt_binary<ET>
{
  Lazy_exact_nt_Min (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_nt_binary<ET>(a,b) {}

  void update_approx() const
  {
      // in = min(op1.approx(), op2.approx());
      in = Min<Interval_nt_advanced>()(op1.approx(), op2.approx());
  }
  ET exact() const
  {
      if (!et)
	  // et = new ET(min(op1.exact(), op2.exact()));
	  et = new ET(Min<ET>()(op1.exact(), op2.exact()));
      return *et;
  }
};

// Maximum
template <typename ET>
struct Lazy_exact_nt_Max : public Lazy_exact_nt_binary<ET>
{
  Lazy_exact_nt_Max (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_nt_binary<ET>(a,b) {}

  void update_approx() const
  {
      in = max(op1.approx(), op2.approx());
  }
  ET exact() const
  {
      if (!et)
	  et = new ET(max(op1.exact(), op2.exact()));
      return *et;
  }
};


// The real number type, handle class
template <typename ET>
struct Lazy_exact_nt : public Handle_for<Lazy_exact_nt_rep<ET> >
{
  typedef Handle_for<Lazy_exact_nt_rep<ET> > PTR;
  typedef Lazy_exact_nt<ET> Self;
  typedef Lazy_exact_nt_rep<ET> Self_rep;
  // typedef Lazy_exact_nt_rep<ET> Self_rep;

  // Data members:
  // Here we should take leda real's idea to cache an exact double value.
  // double d;
  // bool is_double;
  // if "is_double" is false then ptr != NULL.
  // Can we say is_double = (ptr == NULL) ?  And factorize ;-)
  //
  // Computation over this cache can be done cleanly by testing the "inexact"
  // bit of the FPUSW.

  // Ctors:
  Lazy_exact_nt () {}  // Note : this allocates 1 element...

  Lazy_exact_nt (Self_rep *r)
    : PTR(r) {}

  Lazy_exact_nt (const Self & s)
    : PTR(s) {}

#if 0 // provided by Handle_for, normaly...
  Self & operator= (const Self & s)
  {
      Self_rep::inc_count(s.ptr);
      Self_rep::dec_count(ptr);
      ptr = s.ptr;
      return *this;
  }
#endif

  // Operations
  Lazy_exact_nt (const double d)
    : PTR (new Lazy_exact_nt_Cst<ET>(d)) {}
  Lazy_exact_nt (const int i)
    : PTR (new Lazy_exact_nt_Cst<ET>(double(i))) {}

  Self operator+ (const Self & a) const
  { return new Lazy_exact_nt_Add<ET>(*this, a); }

  Self operator- (const Self & a) const
  { return new Lazy_exact_nt_Sub<ET>(*this, a); }

  Self operator* (const Self & a) const
  { return new Lazy_exact_nt_Mul<ET>(*this, a); }

  Self operator/ (const Self & a) const
  { return new Lazy_exact_nt_Div<ET>(*this, a); }

  // Dtor:
  // ~Lazy_exact_nt () {}

  Interval_nt_advanced approx() const
  {
      return ptr->approx();
  }
  ET exact() const
  {
      return ptr->exact();
  }
};

// Other operators are currently provided by STL, the CGAL ones are template.
template <typename ET>
inline
bool
operator< (const Lazy_exact_nt<ET> &r, const Lazy_exact_nt<ET> &s)
{
  // No need for exceptions... could be optimized.
  // Can we have operator< (nothrow), like new (nothrow) ?
  // Could be wonderful...
  try {
    return r.approx() < s.approx();
  }
  catch (...) {
    std::cerr << "Interval filter failure" << std::endl;
    return r.exact() < s.exact();
  }
}

template <typename ET>
inline
Lazy_exact_nt<ET>
abs(const Lazy_exact_nt<ET> & a)
{ return new Lazy_exact_nt_Abs<ET>(a); }

template <typename ET>
inline
Lazy_exact_nt<ET>
sqrt(const Lazy_exact_nt<ET> & a)
{ return new Lazy_exact_nt_Sqrt<ET>(a); }

template <typename ET>
inline
Lazy_exact_nt<ET>
square(const Lazy_exact_nt<ET> & a)
{ return new Lazy_exact_nt_Square<ET>(a); }

template <typename ET>
inline
Lazy_exact_nt<ET>
min(const Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{ return new Lazy_exact_nt_Min<ET>(a, b); }     // Current solution
// { return new Lazy_exact_nt_binary_op<ET,Min>(a, b); } // template² sol.

template <typename ET>
inline
Lazy_exact_nt<ET>
max(const Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{ return new Lazy_exact_nt_Max<ET>(a, b); }

template <typename ET>
std::ostream &
operator<< (std::ostream & os, const Lazy_exact_nt<ET> & I)
{ return os << I.approx(); }

template <typename ET>
inline
Lazy_exact_nt<ET>
operator+=(Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{ return a = a + b; }

template <typename ET>
inline
Lazy_exact_nt<ET>
operator-=(Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{ return a = a - b; }

template <typename ET>
inline
Lazy_exact_nt<ET>
operator*=(Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{ return a = a * b; }

template <typename ET>
inline
Lazy_exact_nt<ET>
operator/=(Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{ return a = a / b; }


// Hackery.
inline
int
convert_from_to (const int&, const Lazy_exact_nt<int> &)
{
    return int();
}

#if !defined(CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION) \
 && !defined(CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION)
struct converter<int, Lazy_exact_nt<int> >
{
    static inline int do_it (const Lazy_exact_nt<int> & z)
    {
        return convert_from_to(int(), z);
    }
};
#endif // CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION


CGAL_END_NAMESPACE

#ifdef CGAL_INTERVAL_ARITHMETIC_H
#include <CGAL/Interval_arithmetic/IA_Lazy_exact_nt.h>
#endif // CGAL_INTERVAL_ARITHMETIC_H

#endif // CGAL_LAZY_EXACT_NT_H
