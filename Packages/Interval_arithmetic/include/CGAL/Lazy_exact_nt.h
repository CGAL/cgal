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
 * Lazy_exact_rep which has pure virtual methods .approx() and .exact().
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
 * - It's not worth optimizing the beast to death, since it'll be obsoleted by
 *   the filtered kernel... one day.
 * - Interval rafinement functionnality ?
 * - The next step will be to replace Cartesian<Lazy_exact_nt<ET> > by a kernel
 *   similar to LOOK in functionality : Lazy_exact_Cartesian<FT,ET>.
 * - Predicates should use the filtered advanced version.
 *   [ done via Filtered_exact<Lazy_exact<X>, X>  :-) ]
 * - Geometric constructions could use the interval_advanced.
 *   [ will be done via the filtered Kernel ]
 * - Separate the handle and the representation in 2 files (?)
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

/*
 * Interface of the rep classes:
 * - .approx()      returns Interval_nt and assumes rounding=nearest.
 *                  [ only called from the handle, and declared in the base ]
 * - .approx_adv()  returns Interval_nt_advanced and assumes rounding=+\inf.
 *                  [ called from .approx(), virtual ]
 * - .exact()       returns ET, if not already done, computes recursively
 *                  virtual
 *
 * - .rafine_approx()  later we can do that (having a birthdate like LOOK ?).
 *
 */

CGAL_BEGIN_NAMESPACE

template <typename ET> class Lazy_exact_nt;

// Abstract base representation class
template <typename ET>
struct Lazy_exact_rep : public Ref_counted
{
  Interval_nt_advanced in;
  ET *et;

  // We do lazy interval computation, so how to mark an uninitialized state ?
  // - Invalid interval [1;0].  Drawback is it's slower to test, and conflicts
  //   with assertions in the generic IA code.
  // - A enum value ?  waste of space...
  // - A particular value (1) for et.
  //   et==NULL => interval, and a fortiori exact, are not computed
  //   et==1    => exact is not computed, but interval is
  Lazy_exact_rep ()
      : in(), et(NULL) {}
      // : in(1,0), et(NULL) {}
  Lazy_exact_rep (const double d)
      : in(d), et((ET*)1) {}

  bool valid_approx() const
  {
    return et!=NULL;
    // return et!=(ET *)1;
    // return is_valid(in);
  }

  void set_valid_approx()
  {
    et = (ET*)1;
  }

  bool valid_exact() const
  {
    return et!=NULL && et!=(ET *)1; // et&~1 != 0
  }

  Interval_nt approx()
  {
      if (!valid_approx())
      {
	  FPU_CW_t backup = FPU_get_and_set_cw(CGAL_FE_UPWARD);
	  update_approx();
	  set_valid_approx();
	  FPU_set_cw(backup);
      }
      return in;
  }

  // update_approx() and approx_adv() must be called with rounding -> +inf.
  // NOTE : those 2 functions seem redundant.  Simplify this.
  Interval_nt_advanced approx_adv()
  {
      // CGAL_assertion(correct rounding)
      if (!valid_approx())
      {
	  update_approx();
	  set_valid_approx();
      }
      return in;
  }

  virtual void update_approx() = 0;
  virtual ET exact() = 0;
  virtual ~Lazy_exact_rep () {};
};

// double constant
template <typename ET>
struct Lazy_exact_Cst : public Lazy_exact_rep<ET>
{
  Lazy_exact_Cst (const double a)
      : Lazy_exact_rep<ET>(a) {}

  void update_approx()
  {
      CGAL_assertion(false);
  }
  ET exact()
  {
      if (!valid_exact())
	  et = new ET(in.inf());
      return *et;
  }
};

// Unary  operations: abs, sqrt, square.
// Binary operations: +, -, *, /, min, max.

// Base unary operation
template <typename ET>
struct Lazy_exact_unary : public Lazy_exact_rep<ET>
{
  const Lazy_exact_nt<ET> op1;

  Lazy_exact_unary (const Lazy_exact_nt<ET> &a)
      : Lazy_exact_rep<ET>(), op1(a) {}
};

// Base binary operation
template <typename ET>
struct Lazy_exact_binary : public Lazy_exact_unary<ET>
{
  const Lazy_exact_nt<ET> op2;

  Lazy_exact_binary (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
      : Lazy_exact_unary<ET>(a), op2(b) {}
};

// Macro for unary operations
#define CGAL_LAZY_UNARY_OP(OP, NAME)                                 \
template <typename ET>                                               \
struct NAME : public Lazy_exact_unary<ET>                            \
{                                                                    \
  NAME (const Lazy_exact_nt<ET> &a)                                  \
      : Lazy_exact_unary<ET>(a) {}                                   \
                                                                     \
  void update_approx()                                               \
  {                                                                  \
      in = OP(op1.approx_adv());                                     \
  }                                                                  \
  ET exact()                                                         \
  {                                                                  \
      if (!valid_exact())                                            \
	  et = new ET(OP(op1.exact()));                              \
      return *et;                                                    \
  }                                                                  \
};

CGAL_LAZY_UNARY_OP(CGAL_NTS abs, Lazy_exact_Abs)
CGAL_LAZY_UNARY_OP(CGAL_NTS square, Lazy_exact_Square)
CGAL_LAZY_UNARY_OP(sqrt, Lazy_exact_Sqrt)

#if 0
// Template binary operator (might be merged with the above ?)
// Note : G++ 2.95 produces an ICE on this  :(  Will try again later...
// Second solution is a macro.
template <typename ET, template <typename T> class Op>
struct Lazy_exact_binary_op : public Lazy_exact_binary<ET>
{
  Lazy_exact_binary_op (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_binary<ET>(a,b) {}

  void update_approx()
  {
      in = Op<Interval_nt_advanced>()(op1.approx_adv(), op2.approx_adv());
  }
  ET exact()
  {
      if (!valid_exact())
	  et = new ET(Op<ET>()(op1.exact(), op2.exact()));
      return *et;
  }
};
#endif

// A macro for +, -, * and /
#define CGAL_LAZY_BINARY_OP(OP, NAME)                                 \
template <typename ET>                                                \
struct NAME : public Lazy_exact_binary<ET>                            \
{                                                                     \
  NAME (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)       \
    : Lazy_exact_binary<ET>(a,b) {}                                   \
                                                                      \
  void update_approx()                                                \
  {                                                                   \
      in = op1.approx_adv() OP op2.approx_adv();                      \
  }                                                                   \
  ET exact()                                                          \
  {                                                                   \
      if (!valid_exact())                                             \
	  et = new ET(op1.exact() OP op2.exact());                    \
      return *et;                                                     \
  }                                                                   \
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
    : Lazy_exact_binary<ET>(a,b) {}

  void update_approx()
  {
      in = min(op1.approx_adv(), op2.approx_adv());
  }
  ET exact()
  {
      if (!valid_exact())
	  et = new ET(min(op1.exact(), op2.exact()));
      return *et;
  }
};

// Maximum
template <typename ET>
struct Lazy_exact_Max : public Lazy_exact_binary<ET>
{
  Lazy_exact_Max (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_binary<ET>(a,b) {}

  void update_approx()
  {
      in = max(op1.approx_adv(), op2.approx_adv());
  }
  ET exact()
  {
      if (!valid_exact())
	  et = new ET(max(op1.exact(), op2.exact()));
      return *et;
  }
};


// The real number type, handle class
template <typename ET>
struct Lazy_exact_nt : public Handle_for<Lazy_exact_rep<ET> >
{
  typedef Handle_for<Lazy_exact_rep<ET> > PTR;
  typedef Lazy_exact_nt<ET> Self;
  typedef Lazy_exact_rep<ET> Self_rep;

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
    : PTR (new Lazy_exact_Cst<ET>(d)) {}
  Lazy_exact_nt (const int i)
    : PTR (new Lazy_exact_Cst<ET>(double(i))) {}

  Self operator+ (const Self & a) const
  { return new Lazy_exact_Add<ET>(*this, a); }

  Self operator- (const Self & a) const
  { return new Lazy_exact_Sub<ET>(*this, a); }

  Self operator* (const Self & a) const
  { return new Lazy_exact_Mul<ET>(*this, a); }

  Self operator/ (const Self & a) const
  { return new Lazy_exact_Div<ET>(*this, a); }

  Interval_nt approx() const
  { return ptr->approx(); }

  Interval_nt_advanced approx_adv() const
  { return ptr->approx_adv(); }

  ET exact() const
  { return ptr->exact(); }

  // The other comparison operators are currently provided by STL.
  bool operator< (const Self & a) const
  {
    // No need for exceptions => could be optimized.
    // Can we have operator< (nothrow), like new (nothrow) ?
    // Could be wonderful...
    try {
      return approx() < a.approx();
    }
    catch (...) {
      std::cerr << "Interval filter failure (<)" << std::endl;
      return exact() < a.exact();
    }
  }
};

namespace NTS {

template <typename ET>
inline
Sign
sign(const Lazy_exact_nt<ET> & a)
{
  try {
    return CGAL_NTS sign(a.approx());
  }
  catch (...) {
    std::cerr << "Interval filter failure (sign)" << std::endl;
    return CGAL_NTS sign(a.exact());
  }
}

template <typename ET>
inline
Comparison_result
compare(const Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{
  try {
    return CGAL_NTS compare(a.approx(), b.approx());
  }
  catch (...) {
    std::cerr << "Interval filter failure (compare)" << std::endl;
    return CGAL_NTS compare(a.exact(), b.exact());
  }
}

template <typename ET>
inline
Lazy_exact_nt<ET>
abs(const Lazy_exact_nt<ET> & a)
{ return new Lazy_exact_Abs<ET>(a); }

template <typename ET>
inline
Lazy_exact_nt<ET>
square(const Lazy_exact_nt<ET> & a)
{ return new Lazy_exact_Square<ET>(a); }

} // namespace NTS

template <typename ET>
inline
Lazy_exact_nt<ET>
sqrt(const Lazy_exact_nt<ET> & a)
{ return new Lazy_exact_Sqrt<ET>(a); }

template <typename ET>
inline
Lazy_exact_nt<ET>
min(const Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{ return new Lazy_exact_Min<ET>(a, b); }
// { return new Lazy_exact_binary_op<ET,Min>(a, b); } // if template² param

template <typename ET>
inline
Lazy_exact_nt<ET>
max(const Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{ return new Lazy_exact_Max<ET>(a, b); }

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

template <typename ET>
inline
io_Operator
io_tag (const Lazy_exact_nt<ET>&)
{ return io_Operator(); }
 
template <typename ET>
inline
Number_tag
number_type_tag (const Lazy_exact_nt<ET>&)
{ return Number_tag(); }

// Temporary hack
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

struct converter<leda_real, Lazy_exact_nt<leda_real> >
{
    static inline leda_real do_it (const Lazy_exact_nt<leda_real> & z)
    {
        return z.exact();
    }
};
#endif // CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION


CGAL_END_NAMESPACE

#ifdef CGAL_INTERVAL_ARITHMETIC_H
#include <CGAL/Interval_arithmetic/IA_Lazy_exact_nt.h>
#endif // CGAL_INTERVAL_ARITHMETIC_H

#endif // CGAL_LAZY_EXACT_NT_H
