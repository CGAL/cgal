// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_LAZY_EXACT_NT_H
#define CGAL_LAZY_EXACT_NT_H

#include <CGAL/basic.h>
#include <CGAL/number_utils.h>
#include <CGAL/Interval_arithmetic.h>

/*
 * This file contains the definition of an interface class: Lazy_exact_nt<ET>.
 * ET is an exact number type (must provide the exact operations needed).
 *
 * Lazy_exact_nt<> provides a DAG-based laxy evaluation, like LEDA's real,
 * Core's Expr, LEA's lazy rationals and a few other ones.
 *
 * The approximation part is based on Interval_nt.
 * The exactness part is provided by ET.
 *
 * The DAG is managed by:
 * - new/delete for the memory part and
 * - virtual functions to denote the various operators.
 */

/*
 * We could also add a string constant to store the variable name, and
 * recursively compute an expression string (a la MetaCGAL).
 */

/*
 * Other packages with similar functionalities:
 * - CGAL::Handle
 * - leda_handle
 * - leda_real
 * - Lazy/LEA
 * - Core/Real/Expr
 * - APU
 * - MetaCGAL
 * - See also C++:
 *   - dynamic types.
 *   - delete order of objects.
 */

/*
 * There are 2 major questions:
 * - How to deal with the DAG management.
 * - How to deal with the dynamic type of a DAG cell.
 */

/*
 * I see only 2 ways to deal with the DAG cell type of operation:
 * - using one enum by operation, and a switch/case.
 * - using one class by operation, and virtual function calls.
 */

/*
 * The decisions should not be based on the predicates, since they will
 * be overloaded to not use any dynamic thing (DAG and memory management).
 * However, it would be nice it we had something that works great in all cases.
 * Would also be nice if the cost was especially low for algorithms that don't
 * do anything outside the predicates (triangulation, convex-hull...).
 */

/*
 * We can also try to manage contructions at the highest possible level.
 * At the NT level, it means all the CGAL utilities (min/max, abs, square...).
 * At the Kernel level: dotproduct ? generic constructions ?
 */

CGAL_BEGIN_NAMESPACE

template <typename ET> class Lazy_exact_nt;

/*
 * The implementation I choose (for the moment) uses virtual functions instead
 * of the enum dirty stuff.
 * The base class only stores a reference counter, and has pure virtual
 * methods:
 * - .interval()
 * - .exact()
 * From this class derives one class per operation, with one constructor.
 *
 * The front class Lazy_exact_nt<> is just a handle to such a class.
 */

/*
 * Other choices that have to be made about .interval():
 * - virtual or not ?
 * - lazy evaluation or not ?
 *
 * lazy        => less rounding mode changes.
 *             => virtual.
 * virtual     => less space for constants (for double we loose only 8bytes)
 *                more flexibility.
 * non virtual => no lazy
 *
 * For the moment, I choose to implement it lazily (=> virtual).
 */

/*
 * Should we have a CT template parameter (with corresponding ctor) ?
 */

/*
 * Do we want to have an interval rafinement functionnality ?
 * (that is: recompute a smaller interval when exact() is called.
 */

/*
 * I see also 2 possibilities for the reference counting:
 * - Use NULL as a pointer to a non valid object, and test for it in a few
 *   places.
 * - Use a static object "lazy_null", with a faked ref counter, to avoid these
 *   tests.  In this case, it's better if it's not template...
 * Benchmark will tell which approach is the best.
 * Note that GCC-2.96 does NULL constant propagation... :)
 */

// unsigned int total_num_objs=0;

// Main base class.
template <typename ET>
class Lazy_exact_nt_dyn_rep
{
public:
  friend Lazy_exact_nt<ET>;
  mutable unsigned int count;
  typedef Lazy_exact_nt_dyn_rep<ET> Self;
public:
  Lazy_exact_nt_dyn_rep () : count(1)
  {
      // total_num_objs++;
      // std::cout << "NEW total num objects = " << total_num_objs << std::endl;
  }
  virtual Interval_nt interval() const = 0;
  virtual ET exact() const = 0;
  virtual ~Lazy_exact_nt_dyn_rep () {};
  static void dec_count(const Self * rep)
  {
    if (rep && --rep->count == 0)
    {
	// total_num_objs--;
	// std::cout << "DELETE total num objects = " << total_num_objs << std::endl;
        delete rep;
    }
  }
  static void inc_count(const Self * rep)
  {
    if (rep)
      rep->count++;
  }
  // virtual void update_interval() const = 0;
  // virtual ostream operator<<() const = 0; // ou string, comme Core ?
};

// double constant.
template <typename ET>
class Lazy_exact_nt_dyn_cst : public Lazy_exact_nt_dyn_rep<ET>
{
public:
  friend Lazy_exact_nt<ET>;
  double d;
  Lazy_exact_nt_dyn_cst (const double a)
      : Lazy_exact_nt_dyn_rep<ET>(), d(a) {}
  Interval_nt interval() const { return d; }
  // void update_interval() const {};
  ET exact() const { return ET(d); }
  ~Lazy_exact_nt_dyn_cst()
  {
      CGAL_assertion(count == 0);
  }
};

// float constant.
template <typename ET>
class Lazy_exact_nt_dyn_float_cst : public Lazy_exact_nt_dyn_rep<ET>
{
public:
  friend Lazy_exact_nt<ET>;
  float d;
  Lazy_exact_nt_dyn_float_cst (const float a)
      : Lazy_exact_nt_dyn_rep<ET>(), d(a) {}
  Interval_nt interval() const { return d; }
  // void update_interval() const {};
  ET exact() const { return ET(d); }
  ~Lazy_exact_nt_dyn_float_cst()
  {
      CGAL_assertion(count == 0);
  }
};

// Unary operations: (probably some factorization of code is welcome...)
// constant, abs, sqrt, square.

// Generic Unary Operation.
template <typename ET>
class Lazy_exact_nt_dyn_unary : public Lazy_exact_nt_dyn_rep<ET>
{
public:
  friend Lazy_exact_nt<ET>;
  mutable ET *et;
  mutable Interval_nt in;
  const Lazy_exact_nt_dyn_rep<ET> *op1;
  Lazy_exact_nt_dyn_unary (const Lazy_exact_nt_dyn_rep<ET> *a)
      : Lazy_exact_nt_dyn_rep<ET>(), et(NULL), in(1,0), op1(a)
  {
      inc_count(a);
  }
  // Interval_nt interval() const {};
  // virtual ET exact() const = 0;
  ~Lazy_exact_nt_dyn_unary()
  {
      CGAL_assertion(count == 0);
      if (et != NULL) // Useless, but faster.
	  delete et;
      dec_count(op1);
  }
};

// abs.
template <typename ET>
class Lazy_exact_nt_dyn_abs : public Lazy_exact_nt_dyn_unary<ET>
{
public:
  friend Lazy_exact_nt<ET>;
  Lazy_exact_nt_dyn_abs (const Lazy_exact_nt_dyn_rep<ET> *a)
      : Lazy_exact_nt_dyn_unary<ET>(a) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = abs(op1->interval());
    return in;
  }
  ET exact() const
  {
      if (et == NULL)
	  et = new ET(abs(op1->exact()));
      return *et;
  }
};

// sqrt.
template <typename ET>
class Lazy_exact_nt_dyn_sqrt : public Lazy_exact_nt_dyn_unary<ET>
{
public:
  friend Lazy_exact_nt<ET>;
  Lazy_exact_nt_dyn_sqrt (const Lazy_exact_nt_dyn_rep<ET> *a)
      : Lazy_exact_nt_dyn_unary<ET>(a) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = sqrt(op1->interval());
    return in;
  }
  ET exact() const
  {
      if (et == NULL)
	  et = new ET(CGAL::sqrt(op1->exact()));
      return *et;
  }
};

// square.
template <typename ET>
class Lazy_exact_nt_dyn_square : public Lazy_exact_nt_dyn_unary<ET>
{
public:
  friend Lazy_exact_nt<ET>;
  Lazy_exact_nt_dyn_square (const Lazy_exact_nt_dyn_rep<ET> *a)
      : Lazy_exact_nt_dyn_unary<ET>(a) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = square(op1->interval());
    return in;
  }
  ET exact() const
  {
      if (et == NULL)
	  et = new ET(square(op1->exact()));
      return *et;
  }
};



// Binary operations: (probably some factorization of code is welcome...)
// +, -, *, /, min, max.

// Generic Binary Operation.  (should it derive from unary instead ?)
template <typename ET>
class Lazy_exact_nt_dyn_binary : public Lazy_exact_nt_dyn_rep<ET>
{
  friend Lazy_exact_nt<ET>;
protected:
  mutable Interval_nt in;
  mutable ET *et;
  const Lazy_exact_nt_dyn_rep<ET> *op1, *op2;
  Lazy_exact_nt_dyn_binary (const Lazy_exact_nt_dyn_rep<ET> *a,
	                    const Lazy_exact_nt_dyn_rep<ET> *b)
      : Lazy_exact_nt_dyn_rep<ET>(), in(1,0), et(NULL), op1(a), op2(b)
  {
      inc_count(a);
      inc_count(b);
  }
  // virtual ET exact() const = 0;
  ~Lazy_exact_nt_dyn_binary()
  {
      CGAL_assertion(count == 0);
      if (et != NULL) // Useless, but faster.
	  delete et;
      dec_count(op1);
      dec_count(op2);
  }
};

// Addition.
template <typename ET>
class Lazy_exact_nt_dyn_add : public Lazy_exact_nt_dyn_binary<ET>
{
  friend Lazy_exact_nt<ET>;
  Lazy_exact_nt_dyn_add (const Lazy_exact_nt_dyn_rep<ET> *a,
	                 const Lazy_exact_nt_dyn_rep<ET> *b)
    : Lazy_exact_nt_dyn_binary<ET>(a,b) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = op1->interval() + op2->interval();
    return in;
  }
  ET exact() const
  {
      if (et == NULL)
	  et = new ET(op1->exact() + op2->exact());
      return *et;
  }
};

// Subtraction.
template <typename ET>
class Lazy_exact_nt_dyn_sub : public Lazy_exact_nt_dyn_binary<ET>
{
  friend Lazy_exact_nt<ET>;
  Lazy_exact_nt_dyn_sub (const Lazy_exact_nt_dyn_rep<ET> *a,
	                 const Lazy_exact_nt_dyn_rep<ET> *b)
    : Lazy_exact_nt_dyn_binary<ET>(a,b) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = op1->interval() - op2->interval();
    return in;
  }
  ET exact() const
  {
      if (et == NULL)
	  et = new ET(op1->exact() - op2->exact());
      return *et;
  }
};

// Multiplication.
template <typename ET>
class Lazy_exact_nt_dyn_mul : public Lazy_exact_nt_dyn_binary<ET>
{
  friend Lazy_exact_nt<ET>;
  Lazy_exact_nt_dyn_mul (const Lazy_exact_nt_dyn_rep<ET> *a,
	                 const Lazy_exact_nt_dyn_rep<ET> *b)
    : Lazy_exact_nt_dyn_binary<ET>(a,b) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = op1->interval() * op2->interval();
    return in;
  }
  ET exact() const
  {
      if (et == NULL)
	  et = new ET(op1->exact() * op2->exact());
      return *et;
  }
};

// Division.
template <typename ET>
class Lazy_exact_nt_dyn_div : public Lazy_exact_nt_dyn_binary<ET>
{
  friend Lazy_exact_nt<ET>;
  Lazy_exact_nt_dyn_div (const Lazy_exact_nt_dyn_rep<ET> *a,
	                 const Lazy_exact_nt_dyn_rep<ET> *b)
    : Lazy_exact_nt_dyn_binary<ET>(a,b) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = op1->interval() / op2->interval();
    return in;
  }
  ET exact() const
  {
      if (et == NULL)
	  et = new ET(op1->exact() / op2->exact());
      return *et;
  }
};

// Minimum.
template <typename ET>
class Lazy_exact_nt_dyn_min : public Lazy_exact_nt_dyn_binary<ET>
{
  friend Lazy_exact_nt<ET>;
  Lazy_exact_nt_dyn_min (const Lazy_exact_nt_dyn_rep<ET> *a,
	                 const Lazy_exact_nt_dyn_rep<ET> *b)
    : Lazy_exact_nt_dyn_binary<ET>(a,b) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = min(op1->interval(), op2->interval());
    return in;
  }
  ET exact() const
  {
      if (et == NULL)
	  et = new ET(min(op1->exact(), op2->exact()));
      return *et;
  }
};

// Maximum.
template <typename ET>
class Lazy_exact_nt_dyn_max : public Lazy_exact_nt_dyn_binary<ET>
{
  friend Lazy_exact_nt<ET>;
  Lazy_exact_nt_dyn_max (const Lazy_exact_nt_dyn_rep<ET> *a,
	                 const Lazy_exact_nt_dyn_rep<ET> *b)
    : Lazy_exact_nt_dyn_binary<ET>(a,b) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = max(op1->interval(), op2->interval());
    return in;
  }
  ET exact() const
  {
      if (et == NULL)
	  et = new ET(max(op1->exact(), op2->exact()));
      return *et;
  }
};



// A few operations are probably still lacking (+=, ... see CGAL doc).
//
// static Lazy_exact_nt_dyn_cst<int> lazy_null(0.0);

// Celui-là devrait contenir un ref_count et tout le bazard sur la pile ?
template <typename ET>
class Lazy_exact_nt
{
public:
  typedef Lazy_exact_nt<ET> Self;
  typedef Lazy_exact_nt_dyn_rep<ET> Self_rep;
  // typedef Lazy_exact_nt_rep<ET> Self_rep;

  // Data member:
  // The rep count could be handled by a non template base class.
  // => non template static object instead of NULL.
  Self_rep *rep;

  Lazy_exact_nt (Self_rep *r)
    : rep(r) {};

public:
  // Ctors:
  Lazy_exact_nt ()
    : rep(NULL) {};
    // : rep(&lazy_null) {}; // Allows to suppress the tests "if (rep)".

  Lazy_exact_nt (const Self & s)
  {
      rep = s.rep;
      Self_rep::inc_count(rep);
  }

  Self & operator= (const Self & s)
  {
      // if (rep != s.rep) // One of these makes the program crash.
      Self_rep::inc_count(s.rep);
      Self_rep::dec_count(rep);
      rep = s.rep;
      return *this;
  }

  // Operations
  Lazy_exact_nt (const double d)
    : rep (new Lazy_exact_nt_dyn_cst<ET>(d)) {}
  Lazy_exact_nt (const int i)
    : rep (new Lazy_exact_nt_dyn_cst<ET>(double(i))) {}
  Lazy_exact_nt (const float i)
    : rep (new Lazy_exact_nt_dyn_float_cst<ET>(i)) {}

  Self operator+ (const Self & a) const
  { return new Lazy_exact_nt_dyn_add<ET>(rep, a.rep); }

  Self operator- (const Self & a) const
  { return new Lazy_exact_nt_dyn_sub<ET>(rep, a.rep); }

  Self operator* (const Self & a) const
  { return new Lazy_exact_nt_dyn_mul<ET>(rep, a.rep); }

  Self operator/ (const Self & a) const
  { return new Lazy_exact_nt_dyn_div<ET>(rep, a.rep); }

  // Dtor:
  ~Lazy_exact_nt ()
  {
      Self_rep::dec_count(rep);
  }

  bool operator< (const Self s) const
  {
    // std::cout << "interval: " << rep->interval() << " <? " << s.rep->interval() << std::endl;
    // std::cout << "exact: " << rep->exact() << " <? " << s.rep->exact() << std::endl;
#if 0
    return rep->interval() < s.rep->interval();
#else
    // No need for exceptions... could be optimized.
    // Can we have operator< (nothrow), like new (nothrow) ?
    // Could be wonderful...
    try {
	return rep->interval() < s.rep->interval();
    }
    catch (...) {
	std::cerr << "Interval filter failure" << std::endl;
	return rep->exact() < s.rep->exact();
    }
#endif
  }
};

template <typename ET>
inline
Lazy_exact_nt<ET>
abs(const Lazy_exact_nt<ET> & a)
{ return new Lazy_exact_nt_dyn_abs<ET>(a.rep); }

template <typename ET>
inline
Lazy_exact_nt<ET>
sqrt(const Lazy_exact_nt<ET> & a)
{ return new Lazy_exact_nt_dyn_sqrt<ET>(a.rep); }

template <typename ET>
inline
Lazy_exact_nt<ET>
square(const Lazy_exact_nt<ET> & a)
{ return new Lazy_exact_nt_dyn_square<ET>(a.rep); }

template <typename ET>
inline
Lazy_exact_nt<ET>
min(const Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{ return new Lazy_exact_nt_dyn_min<ET>(a.rep, b.rep); }

template <typename ET>
inline
Lazy_exact_nt<ET>
max(const Lazy_exact_nt<ET> & a, const Lazy_exact_nt<ET> & b)
{ return new Lazy_exact_nt_dyn_max<ET>(a.rep, b.rep); }

template <typename ET>
std::ostream &
operator<< (std::ostream & os, const Lazy_exact_nt<ET> & I)
{ return os << I.rep->interval(); }


CGAL_END_NAMESPACE

#endif // CGAL_LAZY_EXACT_NT_H
