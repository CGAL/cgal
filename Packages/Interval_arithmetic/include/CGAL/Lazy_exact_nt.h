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
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
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
 * The values are first approximated using Interval_base.
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
 * - Handle_for<> currently doesn't work correctly with class hierarchies.
 *   The de-allocator won't work if it's not the standard one.
 *   Either it is fixed by Stefan, or I must use the usual Handle, or something
 *   else.
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
 * Interface of the rep classes:
 * - .approx()      returns Interval_nt<> (assumes rounding=nearest).
 *                  [ only called from the handle, and declared in the base ]
 * - .exact()       returns ET, if not already done, computes recursively
 *
 * - .rafine_approx()  later we can do that (having a birthdate like LOOK ?).
 *                     could use update_approx().
 */

CGAL_BEGIN_NAMESPACE

template <typename ET> class Lazy_exact_nt;

// Abstract base representation class
template <typename ET>
struct Lazy_exact_rep : public Ref_counted
{
  Interval_base in; // could be const, except for rafinement ? or mutable ?
  ET *et;

  Lazy_exact_rep (const Interval_base i)
      : in(i), et(NULL) {}

  Interval_nt<> approx() const  // Better return a const ref instead ?
  {
      return in;
  }

  ET exact()  // Better return a const ref instead ?
  {
      if (et==NULL)
          update_exact();
      return *et;
  }

  virtual void update_approx() = 0;  // Not used anymore...  at the moment :)
  virtual void update_exact() = 0;
  virtual ~Lazy_exact_rep () {};
};

// double constant
template <typename ET>
struct Lazy_exact_Cst : public Lazy_exact_rep<ET>
{
  Lazy_exact_Cst (const double d)
      : Lazy_exact_rep<ET>(d) {}

  void update_approx() { CGAL_assertion(false); }
  void update_exact()  { et = new ET(in.inf()); }
};

// Unary  operations: abs, sqrt, square.
// Binary operations: +, -, *, /, min, max.

// Base unary operation
template <typename ET>
struct Lazy_exact_unary : public Lazy_exact_rep<ET>
{
  const Lazy_exact_nt<ET> op1;

  Lazy_exact_unary (const Interval_base &i, const Lazy_exact_nt<ET> &a)
      : Lazy_exact_rep<ET>(i), op1(a) {}
};

// Base binary operation
template <typename ET>
struct Lazy_exact_binary : public Lazy_exact_unary<ET>
{
  const Lazy_exact_nt<ET> op2;

  Lazy_exact_binary (const Interval_base &i,
		     const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
      : Lazy_exact_unary<ET>(i, a), op2(b) {}
};

// Here we could use a template class for all operations (STL provides
// function objects plus, minus, multiplies, divides...).  But it would require
// a template template parameter, and GCC 2.95 seems to crash easily with them.

// Macro for unary operations
#define CGAL_LAZY_UNARY_OP(OP, NAME)                                 \
template <typename ET>                                               \
struct NAME : public Lazy_exact_unary<ET>                            \
{                                                                    \
  NAME (const Lazy_exact_nt<ET> &a)                                  \
      : Lazy_exact_unary<ET>(OP(a.approx()), a) {}                   \
                                                                     \
  void update_approx() { in = OP(op1.approx()); }                    \
  void update_exact()  { et = new ET(OP(op1.exact())); }             \
};

CGAL_LAZY_UNARY_OP(CGAL_NTS abs,    Lazy_exact_Abs)
CGAL_LAZY_UNARY_OP(CGAL_NTS square, Lazy_exact_Square)
CGAL_LAZY_UNARY_OP(CGAL::sqrt,      Lazy_exact_Sqrt)

// A macro for +, -, * and /
#define CGAL_LAZY_BINARY_OP(OP, NAME)                                 \
template <typename ET>                                                \
struct NAME : public Lazy_exact_binary<ET>                            \
{                                                                     \
  NAME (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)       \
    : Lazy_exact_binary<ET>(a.approx() OP b.approx(), a, b) {}        \
                                                                      \
  void update_approx() { in = op1.approx() OP op2.approx(); }         \
  void update_exact()  { et = new ET(op1.exact() OP op2.exact()); }   \
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
    : Lazy_exact_binary<ET>(min(a.approx(), b.approx()), a, b) {}

  void update_approx() { in = min(op1.approx(), op2.approx()); }
  void update_exact()  { et = new ET(min(op1.exact(), op2.exact())); }
};

// Maximum
template <typename ET>
struct Lazy_exact_Max : public Lazy_exact_binary<ET>
{
  Lazy_exact_Max (const Lazy_exact_nt<ET> &a, const Lazy_exact_nt<ET> &b)
    : Lazy_exact_binary<ET>(max(a.approx(), b.approx()), a, b) {}

  void update_approx() { in = max(op1.approx(), op2.approx()); }
  void update_exact()  { et = new ET(max(op1.exact(), op2.exact())); }
};


// The real number type, handle class
template <typename ET>
struct Lazy_exact_nt : public Handle_for<Lazy_exact_rep<ET> >  // Buggy :
  // Handle_for<> can't handle virtual hierarchies at the moment.
{
  typedef Handle_for<Lazy_exact_rep<ET> > PTR;
  typedef Lazy_exact_nt<ET> Self;
  typedef Lazy_exact_rep<ET> Self_rep;

  Lazy_exact_nt () {}  // Note : this allocates 1 element

  Lazy_exact_nt (Self_rep *r)
    : PTR(r) {}

  Lazy_exact_nt (const Self & s)
    : PTR(s) {}

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

  Interval_nt<> approx() const  // throw() ?  Can help the compiler...
  { return ptr->approx(); }

  Interval_nt_advanced approx_adv() const
  { return ptr->approx(); }

  ET exact() const
  { return ptr->exact(); }

  // The other comparison operators are currently provided by STL.
  bool operator< (const Self & a) const
  {
    try
    {
      return approx() < a.approx();
    }
    catch (Interval_base::unsafe_comparison)
    {
      std::cerr << "Interval filter failure (<)" << std::endl;
      return exact() < a.exact();
    }
  }
};

// Note:  GCC 2.95 completely and silently ignores the catch block
//        of _template_ function-try-blocks.  GCC 2.96 fixes the bug. 
namespace NTS {

template <typename ET>
inline
Sign
sign(const Lazy_exact_nt<ET> & a)
{
  try
  {
    return CGAL_NTS sign(a.approx());
  }
  catch (Interval_base::unsafe_comparison)
  {
    std::cerr << "Interval filter failure (sign)" << std::endl;
    return CGAL_NTS sign(a.exact());
  }
}

template <typename ET>
inline
Comparison_result
compare(const Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{
  try
  {
    return CGAL_NTS compare(a.approx(), b.approx());
  }
  catch (Interval_base::unsafe_comparison)
  {
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

template <typename ET>
inline
Lazy_exact_nt<ET>
max(const Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{ return new Lazy_exact_Max<ET>(a, b); }

template <typename ET>
std::ostream &
operator<< (std::ostream & os, const Lazy_exact_nt<ET> & a)
{ return os << a.approx(); }

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

#ifndef CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION
template <typename ET>
struct converter<ET, Lazy_exact_nt<ET> >
{
    static inline ET do_it (const Lazy_exact_nt<ET> & z)
    {
        return z.exact();
    }
};
#endif

CGAL_END_NAMESPACE

#ifdef CGAL_INTERVAL_ARITHMETIC_H
#include <CGAL/Interval_arithmetic/IA_Lazy_exact_nt.h>
#endif

#endif // CGAL_LAZY_EXACT_NT_H
