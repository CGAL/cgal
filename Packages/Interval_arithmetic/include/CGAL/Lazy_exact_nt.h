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

// Voir comment font:
// ------------------
// - CGAL::Handle
// - leda_handle
// - leda_real
// - Lazy/LEA
// - Core/Real/Expr
// - APU
// - MetaCGAL
// - C++:
//   - types dynamiques
//   - ordre de deletion des objets

// Deux gros problèmes:
// - mémoire dynamique nécessaire pour le DAG.
// - type dynamique nécessaire pour les opérations.

// Les problèmes arrivent lorsqu'un NT est détruit (ou modifié)
// et que son ref_count est !=0.

// Lazy_exact_nt découpé en 3 morceaux:
// - une partie interne stockant le noeud du DAG
// - une classe allouée dynamiquement.
// - une classe normale, manipulée par le user.

// Deux manières de gérer les operations:
// - enum de l'operation + pointeurs vers les operandes.
// - types dynamiques, avec fonctions virtuelles qui vont bien.

// On se fout un peu du cas des prédicats, ils seront surchargés proprement
// de toute façon.  Ça serait mieux si c'était optimal dans tous les cas, mais
// bon...  Ça serait bien aussi que ça soit optimal à l'exterieur des prédicats
// quand on n'y fait rien...

// Il faut aussi gérer les constructions de plus haut niveau possible:
// min(), max(), abs(), voire d'autres... (dotproduct ?)
// Et les constructions au niveau geometriques ?  Y'a-t-il moyen de les
// integrer a ce niveau facilement ?

CGAL_BEGIN_NAMESPACE

template <typename ET> class Lazy_exact_nt;

/*
template <typename ET>
class Lazy_exact_nt_rep
{
  friend struct Lazy_exact_nt<ET>;
  typedef Lazy_exact_nt_rep<ET> Self;

  // Ctors:
  Lazy_exact_nt_rep () : count(1) {};

  // Status type:
  enum {CST, NEG, ADD, SUB, MUL, DIV} Operation;
  enum {Approx_not_computed=0, Approx_computed=16, Exact_computed=32} Status;

  // Data members:
  Interval	interval;
  Self          *op1, *op2; // Dépend fondamentalement de l'operation...
  ET		*exact;
  int 		status; // = Operation|Status
  unsigned	count;
};
*/

// Do the operators with dynamic typing.
// This class doesn't need to be template (?).
#if 1
template <typename ET>
class Lazy_exact_nt_dyn_rep
{
  friend Lazy_exact_nt<ET>;
  unsigned int count;
public:
  Lazy_exact_nt_dyn_rep () : count(1) {};
  // Interval evaluation:
  // --------------------
  // - virtual function ?
  // - lazy evaluation ?
  // 
  // - lazy        => less rounding changes.
  // - virtual     => less space for constants (for double we loose only 8bytes)
  // - non virtual => no lazy, no virtual calls.
  //
  // (rounding changes are disablable on a whole algorithm basis)
  //  inline ?
  virtual Interval_nt interval() const = 0;
  // virtual ET exact() const = 0;
  // virtual ostream operator<<() const = 0; // ou string, comme Core ?
};
#else
class Lazy_exact_nt_dyn_rep
{
  template <typename ET> Lazy_exact_nt;
  unsigned int count;
  Interval_nt in;
public:
  Lazy_exact_nt_dyn_rep () : count(1), in(1,0) {};
  // Si c'est pas lazy:
  Interval_nt interval() { return in; }
  // Si c'est pas lazy:
  // Interval_nt interval() { if (!in.is_valid()) update_interval(); return in; }
  virtual void update_interval() = 0;
  virtual ET exact() const {}; // = 0;

  // Il faut aussi un destructeur virtuel, pour détruite les *op1, *op2...
  virtual ~Lazy_exact_nt_dyn_rep() {}; // = 0;
};
#endif


// Unary operations: (probably some factorization of code is welcome...)
// constant, abs, sqrt, square.

// double Constant.
template <typename ET>
class Lazy_exact_nt_dyn_cst : public Lazy_exact_nt_dyn_rep<ET>
{
  friend Lazy_exact_nt<ET>;
  const double d;
  // ET *et; // ?
  Lazy_exact_nt_dyn_cst (const double a) : d(a) {}
  Interval_nt interval() const { return d; };
};

// abs.
template <typename ET>
class Lazy_exact_nt_dyn_abs : public Lazy_exact_nt_dyn_rep<ET>
{
  friend Lazy_exact_nt<ET>;
  mutable Interval_nt in;
  const Lazy_exact_nt_dyn_rep<ET> *op;
  // ET *et; // ?
  Lazy_exact_nt_dyn_abs (const Lazy_exact_nt_dyn_rep<ET> *a) : op(a) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = abs(op->interval());
    return in;
  }
};

// sqrt.
template <typename ET>
class Lazy_exact_nt_dyn_sqrt : public Lazy_exact_nt_dyn_rep<ET>
{
  friend Lazy_exact_nt<ET>;
  mutable Interval_nt in;
  const Lazy_exact_nt_dyn_rep<ET> *op;
  // ET *et; // ?
  Lazy_exact_nt_dyn_sqrt (const Lazy_exact_nt_dyn_rep<ET> *a) : op(a) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = sqrt(op->interval());
    return in;
  }
};

// square.
template <typename ET>
class Lazy_exact_nt_dyn_square : public Lazy_exact_nt_dyn_rep<ET>
{
  friend Lazy_exact_nt<ET>;
  mutable Interval_nt in;
  const Lazy_exact_nt_dyn_rep<ET> *op;
  // ET *et; // ?
  Lazy_exact_nt_dyn_square (const Lazy_exact_nt_dyn_rep<ET> *a) : op(a) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = square(op->interval());
    return in;
  }
};


// Binary operations: (probably some factorization of code is welcome...)
// +, -, *, /, min, max.

// Addition.
template <typename ET>
class Lazy_exact_nt_dyn_add : public Lazy_exact_nt_dyn_rep<ET>
{
  friend Lazy_exact_nt<ET>;
  mutable Interval_nt in;
  mutable ET *ex;
  const Lazy_exact_nt_dyn_rep<ET> *op1, *op2;
  // bool interval_computed; // Simulé par (inf > sup), bad = [1;0].
  // bool exact_computed;    // Simulé par (ex == NULL).
  Lazy_exact_nt_dyn_add (const Lazy_exact_nt_dyn_rep<ET> *a,
	                 const Lazy_exact_nt_dyn_rep<ET> *b)
    : in(1,0), ex(NULL), op1(a), op2(b) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = op1->interval() + op2->interval();
    return in;
  }
};

// Subtraction.
template <typename ET>
class Lazy_exact_nt_dyn_sub : public Lazy_exact_nt_dyn_rep<ET>
{
  friend Lazy_exact_nt<ET>;
  mutable Interval_nt in;
  mutable ET *ex;
  const Lazy_exact_nt_dyn_rep<ET> *op1, *op2;
  Lazy_exact_nt_dyn_sub (const Lazy_exact_nt_dyn_rep<ET> *a,
	                 const Lazy_exact_nt_dyn_rep<ET> *b)
    : in(1,0), ex(NULL), op1(a), op2(b) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = op1->interval() - op2->interval();
    return in;
  }
};

// Multiplication.
template <typename ET>
class Lazy_exact_nt_dyn_mul : public Lazy_exact_nt_dyn_rep<ET>
{
  friend Lazy_exact_nt<ET>;
  mutable Interval_nt in;
  mutable ET *ex;
  const Lazy_exact_nt_dyn_rep<ET> *op1, *op2;
  Lazy_exact_nt_dyn_mul (const Lazy_exact_nt_dyn_rep<ET> *a,
	                 const Lazy_exact_nt_dyn_rep<ET> *b)
    : in(1,0), ex(NULL), op1(a), op2(b) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = op1->interval() * op2->interval();
    return in;
  }
};

// Division.
template <typename ET>
class Lazy_exact_nt_dyn_div : public Lazy_exact_nt_dyn_rep<ET>
{
  friend Lazy_exact_nt<ET>;
  mutable Interval_nt in;
  mutable ET *ex;
  const Lazy_exact_nt_dyn_rep<ET> *op1, *op2;
  Lazy_exact_nt_dyn_div (const Lazy_exact_nt_dyn_rep<ET> *a,
	                 const Lazy_exact_nt_dyn_rep<ET> *b)
    : in(1,0), ex(NULL), op1(a), op2(b) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = op1->interval() / op2->interval();
    return in;
  }
};

// min().
template <typename ET>
class Lazy_exact_nt_dyn_min : public Lazy_exact_nt_dyn_rep<ET>
{
  friend Lazy_exact_nt<ET>;
  mutable Interval_nt in;
  mutable ET *ex;
  const Lazy_exact_nt_dyn_rep<ET> *op1, *op2;
  Lazy_exact_nt_dyn_min (const Lazy_exact_nt_dyn_rep<ET> *a,
	                 const Lazy_exact_nt_dyn_rep<ET> *b)
    : in(1,0), ex(NULL), op1(a), op2(b) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = min(op1->interval(), op2->interval());
    return in;
  }
};

// max().
template <typename ET>
class Lazy_exact_nt_dyn_max : public Lazy_exact_nt_dyn_rep<ET>
{
  friend Lazy_exact_nt<ET>;
  mutable Interval_nt in;
  mutable ET *ex;
  const Lazy_exact_nt_dyn_rep<ET> *op1, *op2;
  Lazy_exact_nt_dyn_max (const Lazy_exact_nt_dyn_rep<ET> *a,
	                 const Lazy_exact_nt_dyn_rep<ET> *b)
    : in(1,0), ex(NULL), op1(a), op2(b) {}
  Interval_nt interval() const
  {
    if (!is_valid(in))
      in = max(op1->interval(), op2->interval());
    return in;
  }
};

// A few operations are probably still lacking (+=, ... see CGAL doc).
//
// static Lazy_exact_nt_dyn_cst<int> lazy_null(0.0);

// Celui-là devrait contenir un ref_count et tout le bazard sur la pile ?
template <typename ET>
class Lazy_exact_nt
{
  typedef Lazy_exact_nt<ET> Self;
  typedef Lazy_exact_nt_dyn_rep<ET> Self_rep;
  // typedef Lazy_exact_nt_rep<ET> Self_rep;

  // Data member:
  // The rep count could be handled by a non template base class.
  // => non template static object instead of NULL.
  Self_rep *rep;

  Lazy_exact_nt (Self_rep *r)
    : rep(r) {};

  void inc_count() const
  {
    if (rep)
      rep->count++;
  }

  void dec_count() const
  {
    if (rep)
      if (!--rep->count)
        delete rep;
  }

public:
  // Ctors:
  Lazy_exact_nt ()
    : rep(NULL) {};
    // : rep(&lazy_null) {}; // Allows to suppress the tests "if (rep)".

  Lazy_exact_nt (const Self & s)
  {
    if (rep != s.rep)
    {
	rep = s.rep;
	inc_count();
    }
  }

  Self operator= (const Self s)
  {
      if (rep != s.rep)
      {
	  dec_count();
	  rep = s.rep;
	  inc_count();
      }
  }

  // Operations
  Lazy_exact_nt (const double d)
    : rep (new Lazy_exact_nt_dyn_cst<ET>(d)) {}

  Self operator+ (const Self a) const
  {
    return new Lazy_exact_nt_dyn_add<ET>(rep, a.rep);
  }

  Self operator- (const Self a) const
  {
    return new Lazy_exact_nt_dyn_sub<ET>(rep, a.rep);
  }

  Self operator* (const Self a) const
  {
    return new Lazy_exact_nt_dyn_mul<ET>(rep, a.rep);
  }

  Self operator/ (const Self a) const
  {
    return new Lazy_exact_nt_dyn_div<ET>(rep, a.rep);
  }

  // Dtor:
  ~Lazy_exact_nt () { dec_count(); }

  bool operator< (const Self s) const
  {
#if 1
    return rep->interval() < s.rep->interval();
#else
    // No need for exceptions... could be more optimized.
    // Can we have operator< (nothrow), like new (nothrow) ?
    // Could be wonderful...
    try { return rep->interval() < s.rep->interval(); }
    catch pipo() { return rep->exact() < s.rep->exact(); }
#endif
  }
};

template <typename ET>
inline
Lazy_exact_nt<ET>
sqrt(const Lazy_exact_nt<ET> a)
{
    return new Lazy_exact_nt_dyn_sqrt<ET>(a.rep);
}

template <typename ET>
inline
Lazy_exact_nt<ET>
square(const Lazy_exact_nt<ET> a)
{
    return new Lazy_exact_nt_dyn_square<ET>(a.rep);
}

template <typename ET>
inline
Lazy_exact_nt<ET>
abs(const Lazy_exact_nt<ET> a)
{
    return new Lazy_exact_nt_dyn_abs<ET>(a.rep);
}

template <typename ET>
inline
Lazy_exact_nt<ET>
min(const Lazy_exact_nt<ET> a, const Lazy_exact_nt<ET> b)
{
    return new Lazy_exact_nt_dyn_min<ET>(a.rep, b.rep);
}

template <typename ET>
inline
Lazy_exact_nt<ET>
max(const Lazy_exact_nt<ET> a, const Lazy_exact_nt<ET> b)
{
    return new Lazy_exact_nt_dyn_max<ET>(a.rep, b.rep);
}


CGAL_END_NAMESPACE

#endif // CGAL_LAZY_EXACT_NT_H
