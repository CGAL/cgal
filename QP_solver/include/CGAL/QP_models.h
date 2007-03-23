// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Bernd Gaertner <gaertner@inf.ethz.ch>, Kaspar Fischer

#ifndef CGAL_QP_MODELS_H
#define CGAL_QP_MODELS_H

#include <CGAL/basic.h>
#include <CGAL/iterator.h>
#include <CGAL/algorithm.h>
#include <CGAL/QP_solver/basic.h>
#include <CGAL/QP_solver/iterator.h>
#include <vector> 
#include <map>
#include <iomanip>
#include <istream>
#include <sstream>

// this file defines the following models:
// - Quadratic_program_from_iterators
// - Quadratic_program_from_pointers
// - Quadratic_program
// - Nonngative_quadratic_program_from_iterators
// - Nonengative_quadratic_program_from_pointers
// - Nonengative_quadratic_program_
// - Free_quadratic_program_from_iterators
// - Free_quadratic_program_from_pointers
// - Free_quadratic_program
// - Linear_program_from_iterators
// - Linear_program_from_pointers
// - Linear_program
// - Nonngative_linear_program_from_iterators
// - Nonengative_linear_program_from_pointers
// - Nonengative_linear_program
// - Free_linear_program_from_iterators
// - Free_linear_program_from_pointers
// - Free_linear_program
// - Quadratic_program_from_mps
// - Linear_program_from_mps

// for convenience, every model is actually a model of the
// concept QuadraticProgramInterface:
#define QP_MODEL_ITERATOR_TYPES \
  typedef typename Base::A_iterator A_iterator;\
  typedef typename Base::B_iterator B_iterator;\
  typedef typename Base::R_iterator R_iterator;\
  typedef typename Base::FL_iterator FL_iterator;\
  typedef typename Base::L_iterator L_iterator;\
  typedef typename Base::FU_iterator FU_iterator;\
  typedef typename Base::U_iterator U_iterator;\
  typedef typename Base::D_iterator D_iterator;\
  typedef typename Base::C_iterator C_iterator;\
  typedef typename Base::C_entry C_entry

// macros to ease coding of copying models
// invariant: a pointer is either 0, or it points to valid memory
#define QP_MODEL_MAKE_ABRC \
  if (this->n_ > 0) this->a_it = new NT*[this->n_];\
  for (int j=0; j<this->n_; ++j)\
    if (this->m_ > 0)\
      this->a_it[j] = new NT[this->m_];\
    else\
      this->a_it[j] = 0;\
  if (this->m_ > 0) this->b_it = new NT[this->m_];\
  if (this->m_ > 0) this->r_it = new CGAL::Comparison_result[this->m_];\
  if (this->n_ > 0) this->c_it = new NT[this->n_]

#define QP_MODEL_SET_ABRC \
  for (int j=0; j<this->n_; ++j)\
    CGAL::copy_n (*(a+j), this->m_, this->a_it[j]);\
  CGAL::copy_n (b, this->m_, this->b_it);\
  CGAL::copy_n (r, this->m_, this->r_it);\
  CGAL::copy_n (c, this->n_, this->c_it)

#define QP_MODEL_SET_ABRC_FROM_P \
  for (int j=0; j<this->n_; ++j)\
    CGAL::copy_n (*(p.a()+j), this->m_, this->a_it[j]);\
  CGAL::copy_n (p.b(), this->m_, this->b_it);\
  CGAL::copy_n (p.r(), this->m_, this->r_it);\
  CGAL::copy_n (p.c(), this->n_, this->c_it)

#define QP_MODEL_DELETE_ABRC \
  delete[] this->c_it;\
  delete[] this->r_it;\
  delete[] this->b_it;\
  for (int j=0; j<this->n_; ++j) delete[] this->a_it[j];\
  delete[] this->a_it

#define QP_MODEL_MAKE_BOUNDS \
  if (this->n_ > 0) this->fl_it = new bool[this->n_];\
  if (this->n_ > 0) this->l_it = new NT[this->n_];\
  if (this->n_ > 0) this->fu_it = new bool[this->n_];\
  if (this->n_ > 0) this->u_it = new NT[this->n_];

#define QP_MODEL_SET_BOUNDS \
  CGAL::copy_n (fl, this->n_, this->fl_it);\
  CGAL::copy_n (l, this->n_, this->l_it);\
  CGAL::copy_n (fu, this->n_, this->fu_it);\
  CGAL::copy_n (u, this->n_, this->u_it)

#define QP_MODEL_SET_BOUNDS_FROM_P \
  CGAL::copy_n (p.fl(), this->n_, this->fl_it);\
  CGAL::copy_n (p.l(), this->n_, this->l_it);\
  CGAL::copy_n (p.fu(), this->n_, this->fu_it);\
  CGAL::copy_n (p.u(), this->n_, this->u_it)

#define QP_MODEL_DELETE_BOUNDS \
  delete[] this->u_it;\
  delete[] this->fu_it;\
  delete[] this->l_it;\
  delete[] this->fl_it

#define QP_MODEL_MAKE_D \
  if (this->n_ > 0) this->d_it = new NT*[this->n_];\
  for (int j=0; j<this->n_; ++j) this->d_it[j] = new NT[j+1];

#define QP_MODEL_SET_D \
  for (int j=0; j<this->n_; ++j) CGAL::copy_n (*(d+j), j+1, this->d_it[j])

#define QP_MODEL_SET_D_FROM_P \
  for (int j=0; j<this->n_; ++j) CGAL::copy_n (*(p.d()+j), j+1, this->d_it[j])

#define QP_MODEL_DELETE_D \
  for (int j=0; j<this->n_; ++j) delete[] this->d_it[j];\
  delete[] this->d_it


CGAL_BEGIN_NAMESPACE

// default iterator types to used to make linear / nonnegative models
// conform to QuadraticProgramInterface
template <class Iterator>
class QP_model_default_iterators {
private:
  typedef typename std::iterator_traits<Iterator>::value_type value_type;
public:
  typedef Const_oneset_iterator<value_type> 
  It_1d; // 1-dimensional random access iterator for a constant value
  typedef Const_oneset_iterator<Const_oneset_iterator<value_type> > 
  It_2d; // 2-dimensional random access iterator for a constant value
};

// Quadratic_program_from_iterators
// --------------------------------
// this is the base class for all non-mps models
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename FL_it,  // for finiteness of lower bounds (value type bool)
  typename L_it,   // for lower bounds
  typename FU_it,  // for finiteness of upper bounds (value type bool)
  typename U_it,   // for upper bounds
  typename D_it,   // for quadratic matrix D (rowwise)
  typename C_it >  // for objective function c
class Quadratic_program_from_iterators 
{
public:
  // types
  typedef A_it   A_iterator;
  typedef B_it   B_iterator; 
  typedef R_it   R_iterator;
  typedef FL_it  FL_iterator;
  typedef L_it   L_iterator;
  typedef FU_it  FU_iterator;
  typedef U_it   U_iterator;  
  typedef D_it   D_iterator;
  typedef C_it   C_iterator;
  typedef typename std::iterator_traits<C_it>::value_type C_entry;
protected:
  // data
  int n_;
  int m_;
  A_iterator a_it;
  B_iterator b_it; 
  R_iterator r_it;
  FL_iterator fl_it;
  L_iterator l_it;
  FU_iterator fu_it;
  U_iterator u_it; 
  D_iterator d_it;
  C_iterator c_it;
  C_entry c_0; // constant term
public:
  // construction
  Quadratic_program_from_iterators (
     int n, int m, // number of variables / constraints
     const A_iterator& a, 
     const B_iterator& b,
     const R_iterator& r,
     const FL_iterator& fl,
     const L_iterator& l,
     const FU_iterator& fu,
     const U_iterator& u,
     const D_iterator& d,
     const C_iterator& c,
     const C_entry& c0 = C_entry(0))
    : n_ (n), m_ (m), a_it (a), b_it (b), r_it (r), fl_it (fl), l_it (l), 
      fu_it (fu), u_it (u), d_it (d), c_it (c), c_0 (c0)    
  {}

  // access
  int n() const {return n_;}
  int m() const {return m_;}
  A_iterator a() const {return a_it;}
  B_iterator b() const {return b_it;}
  R_iterator r() const {return r_it;}  
  FL_iterator fl() const {return fl_it;}
  L_iterator l() const {return l_it;}
  FU_iterator fu() const {return fu_it;}
  U_iterator u() const {return u_it;}
  D_iterator d() const {return d_it;}
  C_iterator c() const {return c_it;}
  C_entry c0() const {return c_0;}
};

// corresponding global function make_quadratic_program_from_iterators
// -------------------------------------------------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename FL_it,  // for finiteness of lower bounds (value type bool)
  typename L_it,   // for lower bounds
  typename FU_it,  // for finiteness of upper bounds (value type bool)
  typename U_it,   // for upper bounds
  typename D_it,   // for quadratic matrix D (rowwise)
  typename C_it >  // for objective function c
Quadratic_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, D_it, C_it>
make_quadratic_program_from_iterators (
   int n, int m, 
   const A_it& a, 
   const B_it& b, 
   const R_it& r, 
   const FL_it& fl, 
   const L_it& l,
   const FU_it& fu, 
   const U_it& u, 
   const D_it& d, 
   const C_it& c, 
   typename std::iterator_traits<C_it>::value_type c0 = 
   typename std::iterator_traits<C_it>::value_type(0))
{
  return Quadratic_program_from_iterators
    <A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, D_it, C_it>
    (n, m, a, b, r, fl, l, fu, u, d, c, c0);
}	

// Quadratic_program_from_pointers
// -------------------------------
template <typename NT>
class Quadratic_program_from_pointers : 
  public Quadratic_program_from_iterators 
<NT**, NT*, CGAL::Comparison_result*, 
 bool*, NT*, bool*, NT*, NT**, NT*>
{
private:
  typedef Quadratic_program_from_iterators 
  <NT**, NT*, CGAL::Comparison_result*, 
   bool*, NT*, bool*, NT*, NT**, NT*> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  Quadratic_program_from_pointers 
  (
   int n, int m, // number of variables / constraints
   const A_iterator& a, 
   const B_iterator& b,
   const R_iterator& r,
   const FL_iterator& fl,
   const L_iterator& l,
   const FU_iterator& fu,
   const U_iterator& u,
   const D_iterator& d,
   const C_iterator& c,
   const C_entry& c0 = C_entry(0))
    : Base (n, m, a, b, r, fl, l, fu, u, d, c, c0)
  {}  
};

// Quadratic_program (copies the data)
// -----------------------------------
template <typename NT>
class Quadratic_program :
  public Quadratic_program_from_pointers<NT>
{
private:
  typedef Quadratic_program_from_pointers<NT> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  template <typename A_it, typename B_it, typename R_it, typename FL_it, 
	    typename L_it, typename FU_it, typename U_it, typename D_it, 
	    typename C_it>
  Quadratic_program 
  (
   int n, int m, // number of variables / constraints
   const A_it& a, 
   const B_it& b,
   const R_it& r,
   const FL_it& fl,
   const L_it& l,
   const FU_it& fu,
   const U_it& u,
   const D_it& d,
   const C_it& c,
   const C_entry& c0 = C_entry(0))
    : Base (n, m, 0, 0, 0, 0, 0, 0, 0, 0, 0, c0)
  {
    QP_MODEL_MAKE_ABRC;
    QP_MODEL_SET_ABRC;
    QP_MODEL_MAKE_BOUNDS;
    QP_MODEL_SET_BOUNDS;
    QP_MODEL_MAKE_D;
    QP_MODEL_SET_D;
  }

  Quadratic_program ()
    : Base (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  {}

  Quadratic_program (const Quadratic_program& p)
    : Base (p.n(), p.m(), 0, 0, 0, 0, 0, 0, 0, 0, 0, p.c0())
  {
    QP_MODEL_MAKE_ABRC;
    QP_MODEL_SET_ABRC_FROM_P;
    QP_MODEL_MAKE_BOUNDS;
    QP_MODEL_SET_BOUNDS_FROM_P;
    QP_MODEL_MAKE_D;
    QP_MODEL_SET_D_FROM_P;
  }

  Quadratic_program& operator= (const Quadratic_program& p)
  {
    if (this != &p) {
      // cleanup
      QP_MODEL_DELETE_D;
      QP_MODEL_DELETE_BOUNDS;
      QP_MODEL_DELETE_ABRC;
      // reset
      this->n_ = p.n();
      this->m_ = p.m();
      this->c_0 = p.c0();
      QP_MODEL_MAKE_ABRC;
      QP_MODEL_SET_ABRC_FROM_P;
      QP_MODEL_MAKE_BOUNDS;
      QP_MODEL_SET_BOUNDS_FROM_P;
      QP_MODEL_MAKE_D;
      QP_MODEL_SET_D_FROM_P;
    }
    return *this;
  }

  ~Quadratic_program () 
  {
    QP_MODEL_DELETE_D;
    QP_MODEL_DELETE_BOUNDS;
    QP_MODEL_DELETE_ABRC;
  }
};

// Linear_program_from_iterators
// -----------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename FL_it,  // for finiteness of lower bounds (value type bool)
  typename L_it,   // for lower bounds
  typename FU_it,  // for finiteness of upper bounds (value type bool)
  typename U_it,   // for upper bounds
  typename C_it >  // for objective function c
class Linear_program_from_iterators : 
  public Quadratic_program_from_iterators 
<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it,
 typename QP_model_default_iterators<C_it>::It_2d, C_it>
{
private:
  typedef Quadratic_program_from_iterators 
  <A_it, B_it, R_it, FL_it, L_it, FU_it, U_it,
   typename QP_model_default_iterators<C_it>::It_2d, C_it> Base;
  typedef typename QP_model_default_iterators<C_it>::It_2d Const_D_iterator;
public:
   QP_MODEL_ITERATOR_TYPES;
   Linear_program_from_iterators (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,
		     const FL_iterator& fl,
		     const L_iterator& l,
		     const FU_iterator& fu,
		     const U_iterator& u,
		     const C_iterator& c,
		     const C_entry& c0 = C_entry(0)
		     )
    : Base (n, m, a, b, r, fl, l, fu, u, 
	    Const_D_iterator(C_entry(0)), c, c0)
  {}  
}; 

// corresponding global function make_linear_program_from_iterators
// ----------------------------------------------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename FL_it,  // for finiteness of lower bounds (value type bool)
  typename L_it,   // for lower bounds
  typename FU_it,  // for finiteness of upper bounds (value type bool)
  typename U_it,   // for upper bounds
  typename C_it >  // for objective function c
Linear_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, C_it>
make_linear_program_from_iterators (
   int n, int m, 
   const A_it& a, 
   const B_it& b, 
   const R_it& r, 
   const FL_it& fl, 
   const L_it& l,
   const FU_it& fu, 
   const U_it& u, 
   const C_it& c, 
   typename std::iterator_traits<C_it>::value_type c0 =
   typename std::iterator_traits<C_it>::value_type(0))
{
  return Linear_program_from_iterators
    <A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, C_it>
    (n, m, a, b, r, fl, l, fu, u, c, c0);
}	

// Linear_program_from_pointers
// ----------------------------
template <typename NT>
class Linear_program_from_pointers : 
  public Linear_program_from_iterators 
<NT**, NT*, CGAL::Comparison_result*, bool*, NT*, bool*, NT*, NT*>
{
private:
  typedef  Linear_program_from_iterators 
  <NT**, NT*, CGAL::Comparison_result*, bool*, NT*, bool*, NT*, NT*> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  Linear_program_from_pointers (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,
		     const FL_iterator& fl,
		     const L_iterator& l,
		     const FU_iterator& fu,
		     const U_iterator& u,
		     const C_iterator& c,
		     const C_entry& c0 = C_entry(0)
		     )
    : Base (n, m, a, b, r, fl, l, fu, u, c, c0)
  {}  
};

// Linear_program (copies the data)
// --------------------------------
template <typename NT>
class Linear_program :
  public Linear_program_from_pointers<NT>
{
private:
  typedef Linear_program_from_pointers<NT> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  template <typename A_it, typename B_it, typename R_it, typename FL_it, 
	    typename L_it, typename FU_it, typename U_it, typename C_it>
  Linear_program 
  (
   int n, int m, // number of variables / constraints
   const A_it& a, 
   const B_it& b,
   const R_it& r,
   const FL_it& fl,
   const L_it& l,
   const FU_it& fu,
   const U_it& u,
   const C_it& c,
   const C_entry& c0 = C_entry(0))
    : Base (n, m, 0, 0, 0, 0, 0, 0, 0, 0, c0)
  {
    QP_MODEL_MAKE_ABRC;
    QP_MODEL_SET_ABRC;
    QP_MODEL_MAKE_BOUNDS;
    QP_MODEL_SET_BOUNDS;
  }

  Linear_program ()
    : Base (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  {}

  Linear_program (const Linear_program& p)
    : Base (p.n(), p.m(), 0, 0, 0, 0, 0, 0, 0, 0, p.c0())
  {
    QP_MODEL_MAKE_ABRC;
    QP_MODEL_SET_ABRC_FROM_P;
    QP_MODEL_MAKE_BOUNDS;
    QP_MODEL_SET_BOUNDS_FROM_P;
  }

  Linear_program& operator= (const Linear_program& p)
  {
    if (this != &p) {
      // cleanup
      QP_MODEL_DELETE_BOUNDS;
      QP_MODEL_DELETE_ABRC;
      // reset
      this->n_ = p.n();
      this->m_ = p.m();
      this->c_0 = p.c0();
      QP_MODEL_MAKE_ABRC;
      QP_MODEL_SET_ABRC_FROM_P;
      QP_MODEL_MAKE_BOUNDS;
      QP_MODEL_SET_BOUNDS_FROM_P;
    }
    return *this;
  }

  ~Linear_program () 
  {
    QP_MODEL_DELETE_BOUNDS;
    QP_MODEL_DELETE_ABRC;
  }
};


// Nonnegative_quadratic_program_from_iterators
// --------------------------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename D_it,   // for quadratic matrix D (rowwise)
  typename C_it >  // for objective function c
class Nonnegative_quadratic_program_from_iterators : 
  public Quadratic_program_from_iterators <A_it, B_it, R_it, 
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     D_it, C_it>
{
private:
  typedef  Quadratic_program_from_iterators <A_it, B_it, R_it, 
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     D_it, C_it> Base;
  typedef typename QP_model_default_iterators<bool*>::It_1d Const_FLU_iterator;
  typedef typename QP_model_default_iterators<C_it>::It_1d Const_LU_iterator;
public:
   QP_MODEL_ITERATOR_TYPES;
   Nonnegative_quadratic_program_from_iterators (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,		     
		     const D_iterator& d,
		     const C_iterator& c,
		     const C_entry& c0 = C_entry(0)
		     )
    : Base (n, m, a, b, r, 
	    Const_FLU_iterator(true), Const_LU_iterator(C_entry(0)), 
	    Const_FLU_iterator(false), Const_LU_iterator(C_entry(0)), 
	    d, c, c0)
  {}  
};

// corresponding global function 
// make_nonnegative_quadratic_program_from_iterators
// -------------------------------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename D_it,   // for quadratic matrix D (rowwise)
  typename C_it >  // for objective function c
Nonnegative_quadratic_program_from_iterators
<A_it, B_it, R_it, D_it, C_it>
make_nonnegative_quadratic_program_from_iterators (
   int n, int m, 
   const A_it& a, 
   const B_it& b, 
   const R_it& r, 
   const D_it& d, 
   const C_it& c, 
   typename std::iterator_traits<C_it>::value_type c0 =
   typename std::iterator_traits<C_it>::value_type(0))
{
  return Nonnegative_quadratic_program_from_iterators
    <A_it, B_it, R_it, D_it, C_it>
    (n, m, a, b, r, d, c, c0);
}	

// Nonnegative_Quadratic_program_from_pointers
// -------------------------------------------
template <typename NT>
class Nonnegative_quadratic_program_from_pointers : 
  public Nonnegative_quadratic_program_from_iterators 
<NT**, NT*, CGAL::Comparison_result*, NT**, NT*>
{
private:
  typedef Nonnegative_quadratic_program_from_iterators 
<NT**, NT*, CGAL::Comparison_result*, NT**, NT*> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  Nonnegative_quadratic_program_from_pointers (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,
		     const D_iterator& d,
		     const C_iterator& c,
		     const C_entry& c0 = C_entry(0)
		     )
    : Base (n, m, a, b, r, d, c, c0)
  {}  
};

// Nonnegative_quadratic_program (copies the data)
// ----------------------------------------------
template <typename NT>
class Nonnegative_quadratic_program :
  public Nonnegative_quadratic_program_from_pointers<NT>
{
private:
  typedef Nonnegative_quadratic_program_from_pointers<NT> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  template <typename A_it, typename B_it, typename R_it, typename D_it, 
	    typename C_it>
  Nonnegative_quadratic_program 
  (
   int n, int m, // number of variables / constraints
   const A_it& a, 
   const B_it& b,
   const R_it& r,
   const D_it& d,
   const C_it& c,
   const C_entry& c0 = C_entry(0))
    : Base (n, m, 0, 0, 0, 0, 0, c0)
  {
    QP_MODEL_MAKE_ABRC;
    QP_MODEL_SET_ABRC;
    QP_MODEL_MAKE_D;
    QP_MODEL_SET_D;
  }

  Nonnegative_quadratic_program ()
    : Base (0, 0, 0, 0, 0, 0, 0)
  {}

  Nonnegative_quadratic_program (const Nonnegative_quadratic_program& p)
    : Base (p.n(), p.m(), 0, 0, 0, 0, 0, p.c0())
  {
    QP_MODEL_MAKE_ABRC;
    QP_MODEL_SET_ABRC_FROM_P;
    QP_MODEL_MAKE_D;
    QP_MODEL_SET_D_FROM_P;
  }

  Nonnegative_quadratic_program& operator= 
  (const Nonnegative_quadratic_program& p)
  {
    if (this != &p) {
      // cleanup
      QP_MODEL_DELETE_D;
      QP_MODEL_DELETE_ABRC;
      // reset
      this->n_ = p.n();
      this->m_ = p.m();
      this->c_0 = p.c0();
      QP_MODEL_MAKE_ABRC;
      QP_MODEL_SET_ABRC_FROM_P;
      QP_MODEL_MAKE_D;
      QP_MODEL_SET_D_FROM_P;
    }
    return *this;
  }

  ~Nonnegative_quadratic_program () 
  {
    QP_MODEL_DELETE_D;
    QP_MODEL_DELETE_ABRC;
  }
};


// Free_quadratic_program_from_iterators
// -------------------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename D_it,   // for quadratic matrix D (rowwise)
  typename C_it >  // for objective function c
class Free_quadratic_program_from_iterators : 
  public Quadratic_program_from_iterators <A_it, B_it, R_it, 
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     D_it, C_it>
{
private:
  typedef  Quadratic_program_from_iterators <A_it, B_it, R_it, 
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     D_it, C_it> Base;
  typedef typename QP_model_default_iterators<bool*>::It_1d Const_FLU_iterator;
  typedef typename QP_model_default_iterators<C_it>::It_1d Const_LU_iterator;
public:
   QP_MODEL_ITERATOR_TYPES;
   Free_quadratic_program_from_iterators (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,		     
		     const D_iterator& d,
		     const C_iterator& c,
		     const C_entry& c0 = C_entry(0)
		     )
    : Base (n, m, a, b, r, 
	    Const_FLU_iterator(false), Const_LU_iterator(C_entry(0)), 
	    Const_FLU_iterator(false), Const_LU_iterator(C_entry(0)), 
	    d, c, c0)
  {}  
};
// corresponding global function 
// make_free_quadratic_program_from_iterators
// ------------------------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename D_it,   // for quadratic matrix D (rowwise)
  typename C_it >  // for objective function c
Free_quadratic_program_from_iterators
<A_it, B_it, R_it, D_it, C_it>
make_free_quadratic_program_from_iterators (
   int n, int m, 
   const A_it& a, 
   const B_it& b, 
   const R_it& r, 
   const D_it& d, 
   const C_it& c, 
   typename std::iterator_traits<C_it>::value_type c0 =
   typename std::iterator_traits<C_it>::value_type(0))
{
  return Free_quadratic_program_from_iterators
    <A_it, B_it, R_it, D_it, C_it>
    (n, m, a, b, r, d, c, c0);
}	

// Free_Quadratic_program_from_pointers
// -------------------------------------------
template <typename NT>
class Free_quadratic_program_from_pointers : 
  public Free_quadratic_program_from_iterators 
<NT**, NT*, CGAL::Comparison_result*, NT**, NT*>
{
private:
  typedef Free_quadratic_program_from_iterators 
<NT**, NT*, CGAL::Comparison_result*, NT**, NT*> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  Free_quadratic_program_from_pointers (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,
		     const D_iterator& d,
		     const C_iterator& c,
		     const C_entry& c0 = C_entry(0)
		     )
    : Base (n, m, a, b, r, d, c, c0)
  {}  
};

// Free_quadratic_program (copies the data)
// ---------------------------------------
template <typename NT>
class Free_quadratic_program :
  public Free_quadratic_program_from_pointers<NT>
{
private:
  typedef Free_quadratic_program_from_pointers<NT> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  template <typename A_it, typename B_it, typename R_it, typename D_it, 
	    typename C_it>
  Free_quadratic_program 
  (
   int n, int m, // number of variables / constraints
   const A_it& a, 
   const B_it& b,
   const R_it& r,
   const D_it& d,
   const C_it& c,
   const C_entry& c0 = C_entry(0))
    : Base (n, m, 0, 0, 0, 0, 0, c0)
  {
    QP_MODEL_MAKE_ABRC;
    QP_MODEL_SET_ABRC;
    QP_MODEL_MAKE_D;
    QP_MODEL_SET_D;
  }

  Free_quadratic_program ()
    : Base (0, 0, 0, 0, 0, 0, 0)
  {}

  Free_quadratic_program (const Free_quadratic_program& p)
    : Base (p.n(), p.m(), 0, 0, 0, 0, 0, p.c0())
  {
    QP_MODEL_MAKE_ABRC;
    QP_MODEL_SET_ABRC_FROM_P;
    QP_MODEL_MAKE_D;
    QP_MODEL_SET_D_FROM_P;
  }

  Free_quadratic_program& operator= (const Free_quadratic_program& p)
  {
    if (this != &p) {
      // cleanup
      QP_MODEL_DELETE_D;
      QP_MODEL_DELETE_ABRC;
      // reset
      this->n_ = p.n();
      this->m_ = p.m();
      this->c_0 = p.c0();
      QP_MODEL_MAKE_ABRC;
      QP_MODEL_SET_ABRC_FROM_P;
      QP_MODEL_MAKE_D;
      QP_MODEL_SET_D_FROM_P;
    }
    return *this;
  }

  ~Free_quadratic_program () 
  {
    QP_MODEL_DELETE_D;
    QP_MODEL_DELETE_ABRC;
  }
};
 
// Nonnegative_linear_program_from_iterators
// -----------------------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename C_it >  // for objective function c
class Nonnegative_linear_program_from_iterators : 
  public Quadratic_program_from_iterators <A_it, B_it, R_it, 
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     typename QP_model_default_iterators<C_it>::It_2d, C_it>
{
private:
  typedef Quadratic_program_from_iterators <A_it, B_it, R_it, 
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     typename QP_model_default_iterators<C_it>::It_2d, C_it> Base;
  typedef typename QP_model_default_iterators<bool*>::It_1d Const_FLU_iterator;
  typedef typename QP_model_default_iterators<C_it>::It_1d Const_LU_iterator;
  typedef typename QP_model_default_iterators<C_it>::It_2d Const_D_iterator;
public:
   QP_MODEL_ITERATOR_TYPES;
   Nonnegative_linear_program_from_iterators (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,		     
		     const C_iterator& c,
		     const C_entry& c0 = C_entry(0)
		     )
    : Base (n, m, a, b, r, 
	    Const_FLU_iterator(true), Const_LU_iterator(C_entry(0)), 
	    Const_FLU_iterator(false), Const_LU_iterator(C_entry(0)), 
	    Const_D_iterator(C_entry(0)), c, c0)
  {}  
}; 

// corresponding global function 
// make_nonnegative_linear_program_from_iterators
// ----------------------------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename C_it >  // for objective function c
Nonnegative_linear_program_from_iterators
<A_it, B_it, R_it, C_it>
make_nonnegative_linear_program_from_iterators (
   int n, int m, 
   const A_it& a, 
   const B_it& b, 
   const R_it& r, 
   const C_it& c, 
   typename std::iterator_traits<C_it>::value_type c0 =
   typename std::iterator_traits<C_it>::value_type(0))
{
  return Nonnegative_linear_program_from_iterators
    <A_it, B_it, R_it, C_it>
    (n, m, a, b, r, c, c0);
}
	
// Nonnegative_linear_program_from_pointers
// ----------------------------------------
template <typename NT>
class Nonnegative_linear_program_from_pointers : 
  public Nonnegative_linear_program_from_iterators 
<NT**, NT*, CGAL::Comparison_result*, NT*>
{
private:
  typedef Nonnegative_linear_program_from_iterators 
<NT**, NT*, CGAL::Comparison_result*, NT*> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  Nonnegative_linear_program_from_pointers (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,
		     const C_iterator& c,
		     const C_entry& c0 = C_entry(0)
		     )
    : Base (n, m, a, b, r, c, c0)
  {}  
};

// Nonnegative_linear_program (copies the data)
// --------------------------------------------
template <typename NT>
class Nonnegative_linear_program :
  public Nonnegative_linear_program_from_pointers<NT>
{
private:
  typedef Nonnegative_linear_program_from_pointers<NT> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  template <typename A_it, typename B_it, typename R_it, typename C_it>
  Nonnegative_linear_program 
  (
   int n, int m, // number of variables / constraints
   const A_it& a, 
   const B_it& b,
   const R_it& r,
   const C_it& c,
   const C_entry& c0 = C_entry(0))
    : Base (n, m, 0, 0, 0, 0, c0)
  {
    QP_MODEL_MAKE_ABRC;
    QP_MODEL_SET_ABRC;
  }

  Nonnegative_linear_program ()
    : Base (0, 0, 0, 0, 0, 0)
  {}

  Nonnegative_linear_program (const Nonnegative_linear_program& p)
    : Base (p.n(), p.m(), 0, 0, 0, 0, p.c0())
  {
    QP_MODEL_MAKE_ABRC;
    QP_MODEL_SET_ABRC_FROM_P;
  }

  Nonnegative_linear_program& operator= (const Nonnegative_linear_program& p)
  {
    if (this != &p) {
      // cleanup
      QP_MODEL_DELETE_ABRC;
      // reset
      this->n_ = p.n();
      this->m_ = p.m();
      this->c_0 = p.c0();
      QP_MODEL_MAKE_ABRC;
      QP_MODEL_SET_ABRC_FROM_P;
    }
    return *this;
  }

  ~Nonnegative_linear_program () 
  {
    QP_MODEL_DELETE_ABRC;
  }
};


// Free_linear_program_from_iterators
// ----------------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename C_it >  // for objective function c
class Free_linear_program_from_iterators : 
  public Quadratic_program_from_iterators <A_it, B_it, R_it, 
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     typename QP_model_default_iterators<C_it>::It_2d, C_it>
{
private:
  typedef Quadratic_program_from_iterators <A_it, B_it, R_it, 
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     typename QP_model_default_iterators<bool*>::It_1d, 
     typename QP_model_default_iterators<C_it>::It_1d,
     typename QP_model_default_iterators<C_it>::It_2d, C_it> Base;
  typedef typename QP_model_default_iterators<bool*>::It_1d Const_FLU_iterator;
  typedef typename QP_model_default_iterators<C_it>::It_1d Const_LU_iterator;
  typedef typename QP_model_default_iterators<C_it>::It_2d Const_D_iterator;
public:
   QP_MODEL_ITERATOR_TYPES;
   Free_linear_program_from_iterators (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,		     
		     const C_iterator& c,
		     const C_entry& c0 = C_entry(0)
		     )
    : Base (n, m, a, b, r, 
	    Const_FLU_iterator(false), Const_LU_iterator(C_entry(0)), 
	    Const_FLU_iterator(false), Const_LU_iterator(C_entry(0)), 
	    Const_D_iterator(C_entry(0)), c, c0)
  {}  
}; 

// corresponding global function 
// make_free_linear_program_from_iterators
// ---------------------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename C_it >  // for objective function c
Nonnegative_linear_program_from_iterators
<A_it, B_it, R_it, C_it>
make_free_linear_program_from_iterators (
   int n, int m, 
   const A_it& a, 
   const B_it& b, 
   const R_it& r, 
   const C_it& c, 
   typename std::iterator_traits<C_it>::value_type c0 =
   typename std::iterator_traits<C_it>::value_type(0))
{
  return Free_linear_program_from_iterators
    <A_it, B_it, R_it, C_it>
    (n, m, a, b, r, c, c0);
}
	
// Free_linear_program_from_pointers
// ---------------------------------
template <typename NT>
class Free_linear_program_from_pointers : 
  public Free_linear_program_from_iterators 
<NT**, NT*, CGAL::Comparison_result*, NT*>
{
private:
  typedef Free_linear_program_from_iterators 
<NT**, NT*, CGAL::Comparison_result*, NT*> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  Free_linear_program_from_pointers (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,
		     const C_iterator& c,
		     const C_entry& c0 = C_entry(0)
		     )
    : Base (n, m, a, b, r, c, c0)
  {}  
};

// Free_linear_program (copies the data)
// -------------------------------------
template <typename NT>
class Free_linear_program :
  public Free_linear_program_from_pointers<NT>
{
private:
  typedef Free_linear_program_from_pointers<NT> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  template <typename A_it, typename B_it, typename R_it, typename C_it>
  Free_linear_program 
  (
   int n, int m, // number of variables / constraints
   const A_it& a, 
   const B_it& b,
   const R_it& r,
   const C_it& c,
   const C_entry& c0 = C_entry(0))
    : Base (n, m, 0, 0, 0, 0, c0)
  {
    QP_MODEL_MAKE_ABRC;
    QP_MODEL_SET_ABRC;
  }

  Free_linear_program ()
    : Base (0, 0, 0, 0, 0, 0)
  {}

  Free_linear_program (const Free_linear_program& p)
    : Base (p.n(), p.m(), 0, 0, 0, 0, p.c0())
  {
    QP_MODEL_MAKE_ABRC;
    QP_MODEL_SET_ABRC_FROM_P;
  }

  Free_linear_program& operator= (const Free_linear_program& p)
  {
    if (this != &p) {
      // cleanup
      QP_MODEL_DELETE_ABRC;
      // reset
      this->n_ = p.n();
      this->m_ = p.m();
      this->c_0 = p.c0();
      QP_MODEL_MAKE_ABRC;
      QP_MODEL_SET_ABRC_FROM_P;
    }
    return *this;
  }

  ~Free_linear_program () 
  {
    QP_MODEL_DELETE_ABRC;
  }
};

// QP_from_mps
// --------------------------
namespace QP_from_mps_detail {

  // functor to get the appropriate begin-iterator of a container 
  template<typename Container, typename Iterator, typename HowToBegin>
  struct Begin
    : public std::unary_function< Container, Iterator >
  {
    typedef Iterator result_type;
    result_type operator () ( const Container& v) const 
    { 
      return HowToBegin()(v); 
    }
  };

  // the general representation template
  template <typename IT, typename Use_sparse_representation>
  struct  Rep_selector {};

  // sparse representation; uses a map and "fake" random access iterator
  template <typename IT>
  struct Rep_selector<IT, Tag_true> 
  {
    typedef std::map<unsigned int, IT> Container;
    typedef CGAL::Fake_random_access_const_iterator<Container> Iterator;
    struct HowToBegin {
     Iterator operator() (const Container& v) const
      { return Iterator(&v, IT(0));}
    };
    typedef Begin<Container, Iterator, HowToBegin> Beginner;
  };

  // dense representation; uses vector and its const-iterator
  template <typename IT>
  struct Rep_selector<IT, Tag_false> 
  {
    typedef std::vector<IT> Container;
    typedef typename Container::const_iterator Iterator;
    struct HowToBegin {
      Iterator operator() (const Container& v) const
      { return v.begin();}
    };
    typedef Begin<Container, Iterator, HowToBegin> Beginner;
  };

  template <typename Matrix_iter, typename Beginner, typename IT,
	    typename Is_linear>
  struct D_selector {};

  template <typename Matrix_iter, typename Beginner, typename IT>
  struct D_selector<Matrix_iter, Beginner, IT, Tag_true> // linear case
  {
    typedef typename QP_model_default_iterators<IT*>::It_2d D_iterator;
  };

  template <typename Matrix_iter, typename Beginner, typename IT>
  struct D_selector<Matrix_iter, Beginner, IT, Tag_false> // quadratic case
  {
    typedef CGAL::Join_input_iterator_1<Matrix_iter, Beginner> 
    D_iterator;
  };

} // QP_from_mps_detail

template<typename IT_,  // The input number type: the numbers in the
			// MPS stream are expected to be of this type
			// and are read using IT_'s operator>>.
			// Warning: IT_ must support EXACT division by
			// 2 (so x/2 must be exactly representable in
			// IT_).  The reason is that the QMATRIX
			// section of the MPS format stores the matrix
			// 2*D and the QP-solver needs D, so we need
			// to divide by 2.  Note: If the MPS stream is
			// known to only contain a DMATRIX section,
			// this warning can be neglected because in a
			// DMATRIX section, D is stored (and not 2*D).
			// (See also method D_format_type().)
	 typename Is_linear_,
	                // this checks that there is no DMATRIX
                        // section, and it avoids allocation of
	                // an (n x n) matrix in the non-sparse 
	                // case
	 typename Sparse_D_=Tag_false,
                        // Use a sparse representation for D.
	 typename Sparse_A_=Tag_false>
                        // Use a sparse representation for A. 
class QP_from_mps {
public:
  typedef IT_ IT;
  typedef Is_linear_ Is_linear;
  typedef Sparse_D_   Sparse_D;
  typedef Sparse_A_   Sparse_A;

  // choose representation for columns of A-matrix
  typedef typename QP_from_mps_detail::Rep_selector
  <IT, Sparse_A>::Container A_col; 
  typedef typename QP_from_mps_detail::Rep_selector
  <IT, Sparse_A>::Iterator A_col_iterator;

  // choose representation for rows of D-matrix
  typedef typename QP_from_mps_detail::Rep_selector
  <IT, Sparse_D>::Container D_row; 
  typedef typename QP_from_mps_detail::Rep_selector
  <IT, Sparse_D>::Iterator D_row_iterator;

  // choose Beginners (maps col/row containers to their begin iterators)
  typedef typename QP_from_mps_detail::Rep_selector
  <IT, Sparse_A>::Beginner A_Beginner; 
  typedef typename QP_from_mps_detail::Rep_selector
  <IT, Sparse_D>::Beginner D_Beginner;

  // matrices A and D itself are vectors of columns resp. rows 
  typedef std::vector<A_col>           A_Matrix; 
  typedef std::vector<D_row>           D_Matrix;

  // containers for the other QP data
  typedef std::vector<IT>                      Vector;    // b, c, l, u
  typedef std::vector<bool>                    F_vector;  // fl, fu 
  typedef std::vector<CGAL::Comparison_result> R_vector;  // r

  // QP_model types
  typedef CGAL::Join_input_iterator_1
  <typename A_Matrix::const_iterator, A_Beginner > A_iterator;
  typedef typename Vector::const_iterator B_iterator;
  typedef typename R_vector::const_iterator R_iterator;
  typedef typename F_vector::const_iterator FL_iterator;
  typedef typename F_vector::const_iterator FU_iterator;
  typedef typename Vector::const_iterator L_iterator;
  typedef typename Vector::const_iterator U_iterator;

  // select D_iterator, depending on Is_linear
  typedef typename QP_from_mps_detail::D_selector
  <typename D_Matrix::const_iterator, D_Beginner, IT, Is_linear >::D_iterator
  D_iterator;
 //  typedef CGAL::Join_input_iterator_1
//   <typename D_Matrix::const_iterator, D_Beginner>
//     D_iterator;
  typedef typename Vector::const_iterator C_iterator;
  typedef typename std::iterator_traits<C_iterator>::value_type C_entry;
private:
  typedef std::pair<std::string,unsigned int> String_int_pair;
  typedef std::map<std::string,unsigned int> Index_map;
  typedef std::map<unsigned int, IT> Map;

private:
  const bool has_linear_tag;
  const int verbosity_;
  std::istream& from;
  std::string error_msg;
  bool is_format_okay_;
  bool is_linear_;
  const bool use_CPLEX_convention;
  const IT it0;

  // actual problem data to be read from MPS file
  A_Matrix A_;
  D_Matrix D_;                
  Vector b_;                  
  R_vector row_types_;        
  F_vector fl_, fu_;          
  Vector l_, u_;                                
  Vector c_;                 
  IT c_0;   
 
  std::string D_section;      // name of the section from which D was read
  
  // cached data:
  bool is_nonnegative_cached, is_nonnegative_;

  // further data gathered from MPS file:
  std::string name;     // from the NAME section
  std::string comment_; // first comment in the file, if any
  std::string obj;      // name of the objective "constraint"
  Index_map row_names;
  Index_map duplicated_row_names; // to handle RANGES section
  Index_map var_names;
  std::vector<std::string> var_by_index; // name of i-th column  
  std::vector<std::string> row_by_index; // name of i-th row

  // variables used in token() (see below):
  bool use_put_back_token;
  std::string put_back_token;

public: 
  // unofficial helpers used in master_mps_to_derivatives.C;
  // these should be considered private for CGAL end-users:
  const A_Matrix& A_matrix() { return A_; }
  const Vector& b_vector() { return b_; }
  const R_vector& r_vector() { return row_types_; }
  void add_entry_in_A(unsigned int j, unsigned int i, const IT& val) {
    add_entry_in_A (j, i, val, Sparse_A());
  }
private: // helpers to deal with sparse/dense representation:

  void add_column (const Tag_true /*sparse_A*/) 
  {  
    // append a new column to A_ (it's empty by default)
    A_.push_back(A_col());
  }

  void add_column (const Tag_false /*sparse_A*/) 
  {
    // append a new empty column to A_
    A_.push_back(A_col(row_names.size(), IT(0)));
  }

  void initialize_D(unsigned int dimension, const Tag_true /*sparse_D*/ ) 
  {
    // generate dimension many rows (empty by default)
    for (unsigned int i=0; i<dimension; ++i) 
    D_.push_back(D_row()); 
  }
  
  void initialize_D(unsigned int dimension, const Tag_false /*sparse_D*/)
  {
    // generate dimension many empty rows; this can be very costly
    // if dimension is large
    for (unsigned int i=0; i<dimension; ++i)
    D_.push_back(D_row(dimension,IT(0)));
  }

  void set_entry_in_A(unsigned int j,unsigned int i,const IT& val,
		      const Tag_true /*sparse_A*/)      
  {
    // set element i in column j to val
    if (val != it0) A_[j][i] = val;
  }

  void set_entry_in_A(unsigned int j,unsigned int i,const IT& val,
		      const Tag_false /*sparse_A*/) 
  {
    // set element i in column j to val
    A_[j][i] = val;
  }  
  
  void add_entry_in_A(unsigned int j, unsigned int i, const IT& val,
		      const Tag_true /*sparse_A*/)
  {
    // append value to column j; new element has index i
    set_entry_in_A (j, i, val, Tag_true()); 
  }

  void add_entry_in_A(unsigned int j, unsigned int i, const IT& val,
		      const Tag_false /*sparse_A*/)
  {
    // append value to column j; new element has index i
    CGAL_qpe_precondition (i == A_[j].size());
    A_[j].push_back (val);
  }

  const IT& get_entry_in_A(unsigned int j,unsigned int i,
		      const Tag_true /*sparse_A*/)
  {
    typename Map::const_iterator it = A_[j].find(i);
    if (it != A_[j].end())
    return it->second;
    else
    return it0;
  }

  const IT&  get_entry_in_A(unsigned int j,unsigned int i,
		      const Tag_false /*sparse_A*/)
  {
    return A_[j][i];
  }

  void set_entry_in_D(unsigned int i, unsigned int j, const IT& val,
		      const Tag_true /*sparse_D*/) 
  {
    if (val != it0) D_[i][j] = val;
  }

  void set_entry_in_D(unsigned int i, unsigned int j, const IT& val,
		      const Tag_false /*sparse_D*/)
  {
    D_[i][j] = val;
  }

  const IT& get_entry_in_D(unsigned int i,unsigned int j,
		      const Tag_true /*sparse_D*/)
  {
    typename Map::const_iterator it = D_[i].find(j);
    if (it != D_[i].end())
    return it->second;
    else
    return it0;
  }

  const IT&  get_entry_in_D(unsigned int i,unsigned int j,
		      const Tag_false /*sparse_D*/)
  {
    return D_[i][j];
  }

private: // private helpers:

  // Replaces the first occurrence of '%' in msg by replacement.
  std::string replace1(const std::string& msg,const std::string& replacement)
  {
    std::string result(msg);
    const std::string::size_type pos = result.find('%');
    CGAL_qpe_assertion(pos < result.size());
    result.replace(pos,1,replacement);
    return result;
  }

  bool err(const char* msg) {
    error_msg = msg;
    return false;
  }

  bool err1(const char* msg,
	    const std::string& parameter1) {
    error_msg = replace1(msg,parameter1);
    return false;
  }

  bool err2(const char* msg,
	    const std::string& parameter1,
	    const std::string& parameter2) {
    error_msg = replace1(replace1(msg,parameter1),parameter2);
    return false;
  }

  bool err3(const char* msg,
	    const std::string& parameter1,
	    const std::string& parameter2,
	    const std::string& parameter3) {
    error_msg = 
    replace1(replace1(replace1(msg,parameter1),parameter2),parameter3);
    return false;
  }

  void warn(const std::string& msg) {
    std::cerr << "MPS parser warning: " << msg << '.' << std::endl;
  }

  void warn1(const std::string& msg,const std::string& parameter1) {
    warn(replace1(msg,parameter1));
  }

private: // parsing routines:

  // (Note: returns true iff a line-break was eaten.)
  bool whitespace()
  {
    // support for put_token_back():
    if (use_put_back_token)
      return false;
    bool lineBreakFound = false;

    char c;
    bool in_comment = false; // true iff we are within a comment
    const bool remember_comment = comment_.size() == 0;
    while (from.get(c))
      if (in_comment) {
	if (c!='\r' && c!='\n') {
	  if (remember_comment)
	    comment_.push_back(c);
	} else
	  in_comment = false;
      } else { // not in comment?
	if (!isspace(c)) {
	  if (c!='$' && c!='*') {
	    from.putback(c);
	    break;
	  }
	  in_comment = true;
	  lineBreakFound = true;
	} else {
	  if (c=='\r' || c=='\n')
	    lineBreakFound = true;
	}
      }
    return lineBreakFound;
  }
    
  std::string token() {
    if (use_put_back_token) {
      use_put_back_token = false;
      return put_back_token;
    }
    whitespace();
    std::string token;
    char c;
    while (from.get(c)) {
      if (isspace(c)) {
	from.putback(c);
	break;
      }
      token.push_back(c);
    }
    return token;
  }

  void put_token_back(const std::string& token) {
    CGAL_qpe_assertion(!use_put_back_token);
    use_put_back_token = true;
    put_back_token = token;
  }

  template<typename NumberType>
  bool number(NumberType& entry) {
    // whitespace(); the following >> should care for this
    from >> entry;
    return from.good();
  }

  bool name_section();
  bool rows_section();
  bool columns_section();
  bool rhs_section();
  bool ranges_section();
  bool bounds_section();
  bool qmatrix_section();

public: // methods:
  // Create a quadratic program instance from a stream.
  //
  // If use_CPLEX_convention is true and you set the upper bound of a
  // variable to exactly zero and do not give a lower bound then the
  // lower bound will be set to zero as well.  If use_CPLEX_convention
  // and you set the upper bound of a variable to exactly zero and do
  // not specify a lower bound then the lower bound will be set to
  // -infinity.
  QP_from_mps(std::istream& in,bool use_CPLEX_convention=true,
		  int verbosity=0);

  // Returns true if and only if the instance has been properly
  // constructed (i.e., if the QP could be loaded from the MPS
  // stream).
  bool is_valid() const;

  // Returns the verbosity level with which the instance was
  // constructed.
  int verbosity() const
  {
    return verbosity_;
  }

  // If is_valid() returns false, this routine returns an error
  // string describing the problem that was encountered.
  //
  // Precondition: !is_valid()
  const std::string& error();

  // Returns the first comment that was read from the MPS stream.
  const std::string& comment();
  
  // Returns the number of variables in the QP.
  //
  // Precondition: is_valid()
  int n() const
  {
    CGAL_qpe_assertion(is_valid());
    return var_names.size();
  }

  // Returns the number of constraints in the QP.
  //
  // Precondition: is_valid()
  int m() const
  {
    CGAL_qpe_assertion(is_valid());
    return b_.size(); // row_names doesn't work as RANGES may duplicate rows
  }

  // Returns the problem's name (as specified in the MPS-file).
  //
  // Precondition: is_valid()
  const std::string& problem_name()
  {
    CGAL_qpe_assertion(is_valid());
    return name;
  }

  // Returns the name (as present in the MPS-file) of the i-th variable.
  // Note: this routine has linear complexity.
  const std::string& name_of_variable(int i)
  {
    CGAL_qpe_assertion(0<=i && i<n());
    return var_by_index[i];
  }

  // Returns the section name of the MPS stream in which the D matrix
  // (if any present) was read.
  //
  // Note: The MPS file was originally defined for LP's only, and so
  // there was no way to store the D matrix of a QP. Several
  // extensions exists, including:
  //
  // - CPLEX extension: here, not D but the matrix 2*D is stored in a
  //   so-called "QMATRIX" section.  When the parser encounters such a
  //   section it has to divide the read entries by 2 in order to be
  //   able to return (see D_iterator) the matrix D.  For this, IT
  //   must support division by 2.
  // - QP-solver extension: we also support reading D from a section
  //   called "DMATRIX" which is defined exactly like CPLEX's QMATRIX
  //   section with the only difference that we do in fact store D and
  //   not 2*D.
  //
  // Precondition: !is_linear
  const std::string& D_format_type()
  {
    CGAL_qpe_assertion(!is_linear());
    return D_section;
  }

  // Returns an iterator over the matrix A (as needs to be
  // passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  A_iterator a() const
  {
    CGAL_qpe_assertion(is_valid());
    return A_iterator(A_.begin(), A_Beginner());
  }

  // Returns an iterator over the vector b (as needs to be
  // passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  B_iterator b() const
  {
    CGAL_qpe_assertion(is_valid());
    return b_.begin();
  }

  // Returns an iterator over the vector c (as needs to be
  // passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  C_iterator c() const
  {
    CGAL_qpe_assertion(is_valid());
    return c_.begin();
  }

  // Returns the constant term of the objective function (as
  // needs to be passed to the constructor of class QP_solver).
  // Precondition: is_valid()
  IT c0() const
  {
    CGAL_qpe_assertion(is_valid());
    return c_0;
  }

  // Returns an iterator over the matrix D (as needs to be
  // passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  // it calls one of the following two helpers to decide between
  // the appropriate iterators
  D_iterator d() const
  {
    return d(Is_linear());
  }

private:
  D_iterator d (Tag_true /*is_linear*/) const {
    return D_iterator(it0); 
  }

  D_iterator d(Tag_false /*is_linear*/) const {
    CGAL_qpe_assertion(is_valid());
    return D_iterator(D_.begin(), D_Beginner());
  }
  
  // checks symmetry; if false, it assigns offending i and j 
  bool is_symmetric(Tag_true, unsigned int&, unsigned int&) const;  
  bool is_symmetric(Tag_false, unsigned int&, unsigned int&) const;

public:

  // Returns an iterator (of value-type Row_type) over the constraint
  // types (as needs to be passed to the constructor of class
  // QP_solver).
  //
  // Precondition: is_valid()
  R_iterator r() const
  {
    CGAL_qpe_assertion(is_valid());
    return row_types_.begin();
  }

  // Returns an iterator of Booleans specifying for each variable
  // whether it has a finite lower bound or not (such an iterator
  // needs to be passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  FL_iterator fl() const {
    CGAL_qpe_assertion(is_valid());
    return fl_.begin();
  }

  // Returns an iterator of Booleans specifying for each variable
  // whether it has a finite upper bound or not (such an iterator
  // needs to be passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  FU_iterator fu() const {
    CGAL_qpe_assertion(is_valid());
    return fu_.begin();
  }

  // Returns an iterator of IT's specifying for each variable what its
  // lower bound is (provided it is finite); such an iterator needs to
  // be passed to the constructor of class QP_solver.
  //
  // Precondition: is_valid()
  U_iterator u() const {
    CGAL_qpe_assertion(is_valid());
    return u_.begin();
  }

  // Returns an iterator of IT's specifying for each variable what its
  // upper bound is (provided it is finite); such an iterator needs to
  // be passed to the constructor of class QP_solver.
  //
  // Precondition: is_valid()
  L_iterator l() const {
    CGAL_qpe_assertion(is_valid());
    return l_.begin();
  }

  // Returns true iff the MPS stream could be successfully parsed and
  // the loaded QP instance has a D matrix that is zero (i.e., iff the
  // QP is actually an LP).
  bool is_linear()
  {
    return is_format_okay_ && is_linear_;
  }

  // Returns true iff the MPS stream could be successfully parsed and
  // the loaded QP instance is in standard form.
  bool is_nonnegative()
  {
    if (!is_format_okay_)
      return false;

    if (!is_nonnegative_cached) {
      for (unsigned int i=0; i<var_names.size(); ++i)
	if (fl_[i] == false || l_[i] != 0 ||
	    fu_[i] == true) {
	  is_nonnegative_ = false;
	  return is_nonnegative_;
	}
      is_nonnegative_ = true;
    }
    return is_nonnegative_;
  }
};

template<typename IT_, typename Is_linear_, 
	 typename Sparse_D_,
	 typename Sparse_A_>
std::ostream& operator<<(std::ostream& o,
			 QP_from_mps<IT_, Is_linear_,
			 Sparse_D_,
			 Sparse_A_>& qp);

// Quadratic_program_from_mps, inherits from QP_from_mps
// -----------------------------------------------------
template<typename IT_,
	 typename Sparse_D_ = Tag_false,
	 typename Sparse_A_ = Tag_false>
class Quadratic_program_from_mps : 
  public QP_from_mps <IT_, Tag_false, Sparse_D_, Sparse_A_>
{
private:
  typedef QP_from_mps <IT_, Tag_false, Sparse_D_, Sparse_A_>
  Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  Quadratic_program_from_mps(std::istream& in,bool use_CPLEX_convention=true,
		  int verbosity=0)
    : Base (in, use_CPLEX_convention, verbosity)
  {}
};


// Linear_program_from_mps, inherits from QP_from_mps
// --------------------------------------------------
template<typename IT_,
	 typename Sparse_D_ = Tag_false,
	 typename Sparse_A_ = Tag_false>
class Linear_program_from_mps : 
  public QP_from_mps <IT_, Tag_true, Sparse_D_, Sparse_A_>
{
private:
  typedef QP_from_mps <IT_, Tag_true, Sparse_D_, Sparse_A_>
  Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  Linear_program_from_mps(std::istream& in,bool use_CPLEX_convention=true,
		  int verbosity=0)
    : Base (in, use_CPLEX_convention, verbosity)
  {}
};


CGAL_END_NAMESPACE

#include <CGAL/QP_solver/QP_models_impl.h>

#endif // CGAL_QP_MODELS_H
