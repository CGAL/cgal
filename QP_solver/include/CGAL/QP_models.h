// Copyright (c) 1997-2007  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Bernd Gaertner <gaertner@inf.ethz.ch>, Kaspar Fischer

#ifndef CGAL_QP_MODELS_H
#define CGAL_QP_MODELS_H

#include <CGAL/license/QP_solver.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/iterator.h>
#include <CGAL/algorithm.h>
#include <CGAL/QP_solver/basic.h>
#include <CGAL/QP_solver/functors.h>
#include <CGAL/IO/io.h>
#include <vector>
#include <map>
#include <iomanip>
#include <istream>
#include <sstream>

#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/boost/iterator/transform_iterator.hpp>

// this file defines the following models:
// - Quadratic_program_from_iterators
// - Quadratic_program
// - Quadratic_program_from_mps
// - Nonngative_quadratic_program_from_iterators
// - Linear_program_from_iterators
// - Nonngative_linear_program_from_iterators

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
namespace CGAL {

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
  int get_n() const {return n_;}
  int get_m() const {return m_;}
  A_iterator get_a() const {return a_it;}
  B_iterator get_b() const {return b_it;}
  R_iterator get_r() const {return r_it;}
  FL_iterator get_fl() const {return fl_it;}
  L_iterator get_l() const {return l_it;}
  FU_iterator get_fu() const {return fu_it;}
  U_iterator get_u() const {return u_it;}
  D_iterator get_d() const {return d_it;}
  C_iterator get_c() const {return c_it;}
  C_entry get_c0() const {return c_0;}
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


namespace QP_model_detail {
  // maps a container to its begin-iterator, as specified by HowToBegin
  template<typename Container, typename Iterator, typename HowToBegin>
  struct Begin
    : public CGAL::cpp98::unary_function< Container, Iterator >
  {
    typedef Iterator result_type;
    result_type operator () ( const Container& v) const
    {
      return HowToBegin()(v);
    }
  };
}

// Quadratic_program
// -----------------
// sparse representation, entries can be set one by one, overriding
// defaults;
template <typename NT_>
class Quadratic_program
{
public:
  typedef NT_ NT;
private:
  // Sparse_vectors
  typedef std::map<std::size_t, NT>
  Sparse_vector;
  typedef std::map<std::size_t, CGAL::Comparison_result>
  Sparse_r_vector;
  typedef std::map<std::size_t, bool>
  Sparse_f_vector;

  // Sparse_matrix
  typedef std::vector<Sparse_vector>
  Sparse_matrix;

  // Sparse_vector_iterators
  //typedef CGAL::Fake_random_access_const_iterator<Sparse_vector>
  typedef boost::transform_iterator<CGAL::Map_with_default<Sparse_vector>,
                                    boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t> >
  Sparse_vector_iterator;
  typedef boost::transform_iterator<CGAL::Map_with_default<Sparse_r_vector>,
            boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t> >
  Sparse_r_vector_iterator;
  typedef boost::transform_iterator<CGAL::Map_with_default<Sparse_f_vector>,
                                    boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t> >
  Sparse_f_vector_iterator;

  // Sparse_matrix_iterator
  struct HowToBegin
  {
    Sparse_vector_iterator operator() (const Sparse_vector& v) const
    { return Sparse_vector_iterator
        (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(0),
         CGAL::Map_with_default<Sparse_vector>(&v, NT(0)));}
  };

  typedef QP_model_detail::Begin
  <Sparse_vector, Sparse_vector_iterator, HowToBegin> Beginner;
  typedef boost::transform_iterator
  <Beginner, typename Sparse_matrix::const_iterator>
  Sparse_matrix_iterator;

  // program data; we maintain the invariant that only the
  // non-default elements are stored; this means that entries
  // may get removed again, and n_, m_ might be larger than
  // absolutely needed; we also maintain the invariants that
  // a_matrix and d_matrix always have n_ elements

  int                                  n_;
  int                                  m_;
  Sparse_matrix                        a_matrix;
  Sparse_vector                        b_vector;
  Sparse_r_vector                      r_vector;
  Sparse_f_vector                      fl_vector;
  Sparse_vector                        l_vector;
  Sparse_f_vector                      fu_vector;
  Sparse_vector                        u_vector;
  Sparse_vector                        c_vector;
  Sparse_matrix                        d_matrix;
  NT                                   c0_;

  // default settings
  CGAL::Comparison_result              default_r;   // from constructor
  bool                                 default_fl;  // from constructor
  NT                                   default_l;   // from constructor
  bool                                 default_fu;  // from constructor
  NT                                   default_u;   // from constructor
protected:
  bool                                 is_valid_;
private:
  std::string                          error_msg;

  // methods
  // enlarges a_matrix, d_matrix to size s, under the
  // precondition that the previous sizes were smaller
  void grow_a_d (int s)
  {
    CGAL_qpe_assertion( a_matrix.size() == d_matrix.size() );
    CGAL_qpe_assertion( a_matrix.size() < static_cast<unsigned int>(s));
    for (int k = static_cast<int>(a_matrix.size()); k < s; ++k) {
      a_matrix.push_back(Sparse_vector());
      d_matrix.push_back(Sparse_vector());
    }
  }

public:
  // interface types
  typedef Sparse_matrix_iterator                     A_iterator;
  typedef Sparse_vector_iterator                     B_iterator;
  typedef Sparse_r_vector_iterator                   R_iterator;
  typedef Sparse_f_vector_iterator                   FL_iterator;
  typedef Sparse_vector_iterator                     L_iterator;
  typedef Sparse_f_vector_iterator                   FU_iterator;
  typedef Sparse_vector_iterator                     U_iterator;
  typedef Sparse_vector_iterator                     C_iterator;
  typedef Sparse_matrix_iterator                     D_iterator;
  typedef NT                                         C_entry;

  // concept methods
  int get_n() const
  {
    CGAL_qpe_assertion(is_valid());
    return n_;
  }
  int get_m() const
  {
    CGAL_qpe_assertion(is_valid());
    return m_;
  }
  A_iterator get_a() const
  {
    CGAL_qpe_assertion(is_valid());
    return A_iterator (a_matrix.begin(), Beginner());
  }
  B_iterator get_b() const
  {
    CGAL_qpe_assertion(is_valid());
    return B_iterator (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(0),
                        CGAL::Map_with_default<Sparse_vector>
                       (&b_vector, NT(0)));
  }
  R_iterator get_r() const
  {
    CGAL_qpe_assertion(is_valid());
    return R_iterator (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(0),
                        CGAL::Map_with_default<Sparse_r_vector>
                       (&r_vector, default_r));
  }
  FL_iterator get_fl() const
  {
    CGAL_qpe_assertion(is_valid());
    return FL_iterator (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(0),
                        CGAL::Map_with_default<Sparse_f_vector>
                        (&fl_vector, default_fl));
  }
  L_iterator get_l() const
  {
    CGAL_qpe_assertion(is_valid());
    return L_iterator (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(0),
                        CGAL::Map_with_default<Sparse_vector>
                       (&l_vector, default_l));
  }
  FU_iterator get_fu() const
  {
    CGAL_qpe_assertion(is_valid());
    return FU_iterator (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(0),
                        CGAL::Map_with_default<Sparse_f_vector>
                        (&fu_vector, default_fu));
  }
  U_iterator get_u() const
  {
    CGAL_qpe_assertion(is_valid());
    return U_iterator (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(0),
                        CGAL::Map_with_default<Sparse_vector>
                       (&u_vector, default_u));
  }
  C_iterator get_c() const
  {
    CGAL_qpe_assertion(is_valid());
    return C_iterator (boost::counting_iterator<std::size_t,boost::use_default,std::ptrdiff_t>(0),
                        CGAL::Map_with_default<Sparse_vector>
                       (&c_vector, NT(0)));
  }
  D_iterator get_d() const
  {
    CGAL_qpe_assertion(is_valid());
    return D_iterator (d_matrix.begin(), Beginner());
  }
  C_entry get_c0() const
  {
    CGAL_qpe_assertion(is_valid());
    return c0_;
  }

  bool is_valid() const
  {
    return is_valid_;
  }

  const std::string& get_error() const
  {
    CGAL_qpe_assertion(!is_valid());
    return error_msg;
  }

  // default constructor
  Quadratic_program
  (CGAL::Comparison_result relation = CGAL::EQUAL,
   bool finite_lower = true,
   NT lower = 0,
   bool finite_upper = false,
   NT upper = 0)
    : n_(0), m_(0), c0_(0),
      default_r(relation), default_fl(finite_lower),
      default_l(lower), default_fu(finite_upper),
      default_u(upper), is_valid_(true)
  {
    CGAL_qpe_assertion(!finite_lower || !finite_upper || lower <= upper);
  }


  // constructor from iterators
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
   const typename std::iterator_traits<C_it>::value_type& c0 = 0)
    : n_(0), m_(0), c0_(0),
      default_r(CGAL::EQUAL), default_fl(true),
      default_l(0), default_fu(false), default_u(0), is_valid_(true)
  {
    // now copy, using the set methods
    for (int j=0; j<n; ++j) {
      for (int i=0; i<m; ++i)
        set_a (j, i, (*(a+j))[i]);
      set_l (j, *(fl+j), *(l+j));
      set_u (j, *(fu+j), *(u+j));
      set_c (j, *(c+j));
    }
    for (int i=0; i<m; ++i) {
      set_b (i, *(b+i));
      set_r (i, *(r+i));
    }
    for (int i=0; i<n; ++i)
      for (int j=0; j<=i; ++j)
        set_d (i, j, (*(d+i))[j]);
    set_c0 (c0);
  }

  // type of problem
  bool is_linear() const
  {
    CGAL_qpe_assertion(d_matrix.size() == (unsigned int)(n_));
    for (int i=0; i<n_; ++i)
      if (!d_matrix[i].empty()) return false;
    return true;
  }

private:
  // helpers to determine bound status
  // checks whether all bounds in flu are as given by the parameter "finite"
  // default_flu is the default-value of the underlying map that is not
  // stored
  bool all_bounds_are
  (bool finite, const Sparse_f_vector& flu, bool default_flu) const
  {
    if (finite == default_flu)
      return flu.empty();
    else
      // are there exactly n non-default values "finite"?
      return flu.size() == (unsigned int)(n_);
  }

  bool all_bounds_are_zero
  // checks whether all bounds in lu are 0. default_lu is the default-value
  // of the underlying map that is not stored
  (const Sparse_vector& lu, const NT& default_lu) const
  {
    if (CGAL::is_zero(default_lu))
      return lu.empty();
    else {
      // are there exactly n non-default values?
      if (lu.size() != (unsigned int)(n_)) return false;
      // ok, we have to test each of the non-default values against zero
      for (typename Sparse_vector::const_iterator
             j = lu.begin(); j != lu.end(); ++j)
        if (!CGAL::is_zero(j->second)) return false;
      return true;
    }
  }

public:

  bool is_nonnegative() const
  {
    return
      all_bounds_are (true, fl_vector, default_fl) &&
      all_bounds_are_zero (l_vector, default_l) &&
      all_bounds_are (false, fu_vector, default_fu);
  }

  bool is_nonpositive() const
  {
     return
      all_bounds_are (false, fl_vector, default_fl) &&
      all_bounds_are_zero (u_vector, default_u) &&
      all_bounds_are (true, fu_vector, default_fu);
  }

  bool is_free() const
  {
    return
      all_bounds_are (false, fl_vector, default_fl) &&
      all_bounds_are (false, fu_vector, default_fu);
  }

  // set methods
  void set_a (int j, int i, const NT& val)
  {
    CGAL_qpe_assertion (j >= 0);
    CGAL_qpe_assertion (i >= 0);
    if (j >= n_) {
      n_ = j+1;
      grow_a_d(n_);
    }
    if (i >= m_) m_ = i+1;
    if (CGAL::is_zero(val))
      a_matrix[j].erase(i);
    else
      a_matrix[j][i] = val;
  }

  void set_b (int i,  const NT& val)
  {
    CGAL_qpe_assertion (i >= 0);
    if (i >= m_) m_ = i+1;
    if (CGAL::is_zero(val))
      b_vector.erase(i);
    else
      b_vector[i] = val;
  }

  void set_r (int i, CGAL::Comparison_result val)
  {
    CGAL_qpe_assertion (i >= 0);
    if (i >= m_) m_ = i+1;
    if (val == default_r)
      r_vector.erase(i);
    else
      r_vector[i] = val;
  }

  void set_l (int j, bool is_finite, const NT& val = NT(0))
  {
    CGAL_qpe_assertion (j >= 0);
    if (j >= n_) {
      n_ = j+1;
      grow_a_d(n_);
    }
    if (is_finite == default_fl)
      fl_vector.erase(j);
    else
      fl_vector[j] = is_finite;
    if (val == default_l)
      l_vector.erase(j);
    else
      l_vector[j] = val;
  }

  void set_u (int j, bool is_finite, const NT& val = NT(0))
  {
    CGAL_qpe_assertion (j >= 0);
    if (j >= n_) {
      n_ = j+1;
      grow_a_d(n_);
    }
      if (is_finite == default_fu)
      fu_vector.erase(j);
    else
      fu_vector[j] = is_finite;
    if (val == default_u)
      u_vector.erase(j);
    else
      u_vector[j] = val;
  }

  void set_c (int j, const NT& val)
  {
    CGAL_qpe_assertion (j >= 0);
    if (j >= n_) {
      n_ = j+1;
      grow_a_d(n_);
    }
    if (CGAL::is_zero(val))
      c_vector.erase(j);
    else
      c_vector[j] = val;
  }

  void set_c0 (const NT& val)
  {
    c0_ = val;
  }

  void set_d (int i, int j, const NT& val)
  {
    CGAL_qpe_assertion (i >= 0);
    CGAL_qpe_assertion (j >= 0);
    CGAL_qpe_assertion (j <= i); // lower-diagonal entry
    if (i >= n_) {
      n_ = i+1;
      grow_a_d(n_);
    }
    if (CGAL::is_zero(val))
      d_matrix[i].erase(j);
    else
      d_matrix[i][j] = val;
  }

protected:
  // helpers for errors/warnings
  std::string replace1
  (const std::string& msg,const std::string& replacement) const
  {
    std::string result(msg);
    const std::string::size_type pos = result.find('%');
    CGAL_qpe_assertion(pos < result.size());
    result.replace(pos,1,replacement);
    return result;
  }

  bool err(const char* msg) {
    error_msg = msg;
    is_valid_ = false;
    return false;
  }

  bool err1(const char* msg,
            const std::string& parameter1) {
    error_msg = replace1(msg,parameter1);
    is_valid_ = false;
    return false;
  }

  bool err2(const char* msg,
            const std::string& parameter1,
            const std::string& parameter2) {
    error_msg = replace1(replace1(msg,parameter1),parameter2);
    is_valid_ = false;
    return false;
  }

  bool err3(const char* msg,
            const std::string& parameter1,
            const std::string& parameter2,
            const std::string& parameter3) {
    error_msg =
    replace1(replace1(replace1(msg,parameter1),parameter2),parameter3);
    is_valid_ = false;
    return false;
  }

  void warn(const std::string& msg) const {
    std::cerr << "Warning: " << msg << '.' << std::endl;
  }

  void warn1(const std::string& msg,const std::string& parameter1) const {
    warn(replace1(msg,parameter1));
  }


};

// Quadratic_program_from_mps
// --------------------------
// for reading a QP from a stream in MPS format

template <typename NT_>
class Quadratic_program_from_mps :
  public Quadratic_program<NT_>
{
public:
  typedef NT_ NT;
private:
  typedef Quadratic_program<NT> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
private:
  // types
  typedef std::map<std::string,unsigned int> Index_map;
  typedef std::pair<std::string,unsigned int> String_int_pair;
public:
  // constructor
  Quadratic_program_from_mps
  (std::istream& in)
    : Base(), from(in), nt0(0), use_put_back_token(false)
  {
    // read NAME section:
    if (!name_section())
      return;

    // read ROWS section:
    if (!rows_section())
      return;

    // read COLUMNS section:
    if (!columns_section())
      return;

    // read RHS section:
    if (!rhs_section())
      return;

    // check for (optional) RANGES section:
    if (!ranges_section())
      return;

    // read optional BOUNDS section:
    if (!bounds_section())
      return;

    // read optional QMATRIX section:
    if (!qmatrix_section())
      return;

    // check for ENDATA:
    const std::string end = token();
    if (end != "ENDATA") {
      this->err1("ENDDATA expected but found '%'",end);
      return;
    }

    // remember the number of variables/constraint that we have now
    n_after_construction = this->get_n();
    m_after_construction = this->get_m();
  }

  // returns the first comment that was read from the MPS stream
  const std::string& get_comment() const
  {
    return comment_;
  }

  // returns name of the problem
  const std::string& get_problem_name() const
  {
    return name;
  }

  const std::string& variable_name_by_index(int j) const
  {
    CGAL_qpe_assertion(this->is_valid());
    CGAL_qpe_assertion(0<=j && j<n_after_construction);
    return var_by_index[j];
  }

  int variable_index_by_name (const std::string& name) const
  {
    const Index_map::const_iterator var_name = var_names.find(name);
    if (var_name == var_names.end()) // unknown variable
      return -1;
    else
      return var_name->second;
  }

  const std::string& constraint_name_by_index(int i) const
  {
    CGAL_qpe_assertion(this->is_valid());
    CGAL_qpe_assertion(0<=i && i<m_after_construction);
    return row_by_index[i];
  }

  int constraint_index_by_name (const std::string& name) const
  {
    const Index_map::const_iterator row_name = row_names.find(name);
    if (row_name == row_names.end()) // unknown constraint
      return -1;
    else
      return row_name->second;
  }

private:
  // data
  // ----
  std::istream& from;    // input stream
  NT nt0;

  std::string D_section; // name of the section from which D was read
  std::string name;      // name of the problem
  std::string comment_;  // first comment in the input, if any
  std::string obj;       // name of the objective "constraint"
  int n_after_construction;
  int m_after_construction;

  Index_map row_names;
  Index_map duplicated_row_names; // to handle RANGES section
  Index_map var_names;
  std::vector<std::string> var_by_index; // name of i-th column
  std::vector<std::string> row_by_index; // name of i-th row

  // variables used in token() (see below):
  bool use_put_back_token;
  std::string put_back_token;

  // parsing routines
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
    from >> CGAL::iformat(entry);
    return from.good();
  }

  bool name_section()
  {
    const std::string t = token();
    if (t != "NAME")
      return this->err("expected 'NAME'");
    // NAME: everything found until line break; whitespaces are allowed
    char c;
    std::string token;
    std::string whitespaces;
    if (whitespace())
      // line break eaten, name is empty
      return true;
    do {
      from.get(c);
      if (c == '\r' || c == '\n') break;
      if (isspace(c))
        whitespaces.push_back(c); // save whitespace
      else {
        // new actual character found: previous whitespaces belong to name
        name += whitespaces;
        whitespaces.clear();
        name.push_back(c);
      }
    } while (true);
    return true;
  }

  bool rows_section()
  {
    std::string t = token();
    if (t != "ROWS")
      return this->err1("expected 'ROWS' but found '%'",t);

    // read 'N', 'G', 'L', or 'E', and the name of the constraint:
    t = token();
    int i = 0; // row index
    while (t != "COLUMNS") {
      const char type = t[0];
      const std::string symbol(t); // for error message below
      t = token();
      switch (type) {
      case 'N':
        // register name of objective row:
        if (obj.size() == 0) // remember first (and ignore others)
          obj = t;
        break;
      case 'G':
      case 'L':
      case 'E':
        {
          // register name of >=, <=, or = constraint:
          this->set_r (i,
                 type == 'G'? CGAL::LARGER :
                 (type == 'E'? CGAL::EQUAL : CGAL::SMALLER));
          if (row_names.find(t) != row_names.end())
            return this->err1("duplicate row name '%' in section ROWS",t);
          row_names.insert(String_int_pair(t,i));
          row_by_index.push_back(t);
          ++i;
        }
        break;
      default:
        return
          this->err1
          ("expected 'N', 'L', 'E', or 'G' in ROWS section but found '%'",
           symbol);
      }
      t = token();
    }
    put_token_back(t);

    return true;
  }

  bool columns_section()
  {
    std::string t = token();
    if (t != "COLUMNS")
      return this->err1("expected 'COLUMNS' but found '%'",t);

    t = token();
    while (t != "RHS") {
      // find variable name:
      unsigned int var_index;
      std::string col_name;
      const Index_map::const_iterator var_name = var_names.find(t);
      if (var_name == var_names.end()) { // new variable?
        var_index = static_cast<unsigned int>(var_names.size());
        col_name = t;
        var_names.insert(String_int_pair(t,var_index));
        var_by_index.push_back(t);
      } else { // variable that is already known?
        var_index = var_name->second;
        col_name = var_name->first;
      }

      bool doItAgain = true;
      for (int i=0; doItAgain; ++i) {
        // read row identifier:
        t = token();

        // read number:
        NT val;
        if (!number(val))
          return this->err2
            ("number expected after row identifier '%' in '%' COLUMNS record",
             t,col_name);

        // store number:
        if (t == obj) { // objective row?
          this->set_c(var_index, val);
        } else { // not objective row?
          const Index_map::const_iterator row_name = row_names.find(t);
          if (row_name == row_names.end())
            return this->err1
              ("unknown row identifier '%' in section COLUMNS",t);
          this->set_a (var_index, row_name->second, val);
        }

        // determine if we need to read another number:
        doItAgain = i==0 && !whitespace();
      }

      // read next token:
      t = token();
    }
    put_token_back(t);

    return true;
  }

  bool rhs_section()
  {
    this->set_c0 (nt0);  // no constant term yet
    std::string t = token();
    if (t != "RHS")
      return this->err1("expected 'RHS' but found '%'",t);

    t = token();
    std::string rhs_id;
    while (t != "RANGES" && t != "BOUNDS" &&
           t != "DMATRIX" && t != "QMATRIX" && t != "QUADOBJ" &&
           t != "ENDATA") {
      // read rhs identifier and if it is different from the one
      // from the previous iteration, ignore the whole row:
      bool ignore = false;
      std::string ignored;
      if (rhs_id.size() == 0) { // first time we enter the loop?
        rhs_id = t;
      } else {                  // rhs_id already set
        if (t != rhs_id) {
          ignore = true;        // ignore other rhs identifiers
          ignored = t;
        }
      }

      bool doItAgain = true;
      for (int i=0; doItAgain; ++i) {
        // read row identifier:
        t = token();

        // read number:
        NT val;
        if (!number(val))
          return this->err1("number expected after '%' in this RHS record",t);

        // store number:
        const Index_map::const_iterator row_name = row_names.find(t);
        if (row_name == row_names.end()) {
          // no corresponding constraint; is it the constant term?
          if (t == obj)
            this->set_c0(-val);
          else
            return this->err1("unknown row identifier '%' in section RHS",t);
        } else {
          // we have an actual constraint
          if (!ignore) {
            this->set_b(row_name->second, val);
          } else {
            this->warn1("rhs with identifier '%' ignored", ignored);
          }
        }

        // determine if we need to read another number:
        doItAgain = i==0 && !whitespace();
      }

      // read next token:
      t = token();
    }
    put_token_back(t);

    return true;
  }

  bool ranges_section()
  {
    std::string t = token();
    if (t != "RANGES") { // (Note: RANGES section is optional.)
      put_token_back(t);
      return true;
    }

    t = token();
    std::string range_id;
    while ((t != "BOUNDS" && t != "QMATRIX" &&
            t != "DMATRIX" && t != "QUADOBJ" && t != "ENDATA")) {
      // read rhs identifier and if it is different from the one
      // from the previous iteration, ignore the whole row:
      bool ignore = false;
      std::string ignored;
      if (range_id.size() == 0) { // first time we enter the loop?
        range_id = t;
      } else {                    // range_id already set
        if (t != range_id) {
          ignore = true;          // ignore other range identifiers
          ignored = t;
        }
      }
      bool doItAgain = true;
      for (int i=0; doItAgain; ++i) {
        // read row identifier:
        t = token();

        // read number:
        NT val;
        if (!number(val))
          return this->err1
            ("number expected after '%' in this RANGES record",t);

        // duplicate the constraint, depending on sign of val and type
        // of constraint
        const Index_map::const_iterator row_name = row_names.find(t);
        if (row_name == row_names.end()) {
          return this->err1("unknown row identifier '%' in section RANGES",t);
        } else {
          if (!ignore) {
            int index = row_name->second;
            CGAL::Comparison_result type = *(this->get_r()+index);
            // duplicate the row, unless it has already been duplicated
            const Index_map::const_iterator duplicated_row_name =
              duplicated_row_names.find(t);
            if (duplicated_row_name != duplicated_row_names.end())
              return this->err1
                ("duplicate row identifier '%' in section RANGES",t);
            duplicated_row_names.insert(*row_name);
            std::string dup_name = row_name->first+std::string("_DUPLICATED");
            int new_index = this->get_m();
            row_names.insert(String_int_pair (dup_name, new_index));
            row_by_index.push_back (dup_name);
            for (unsigned int j=0; j<var_names.size(); ++j) {
              NT val = (*(this->get_a()+j))[index];
              this->set_a (j, new_index, val);
            }
            // determine rhs for this new row. Here are the rules:
            // if r is the ranges value and b is the old right-hand
            // side, then we have h <= constraint <= u according to
            // this table:
            //
            // row type       sign of r       h          u
            // ----------------------------------------------
            // G            + or -         b        b + |r|
            // L            + or -       b - |r|      b
            // E              +            b        b + |r|
            // E              -          b - |r|      b

            switch (type) {
            case CGAL::LARGER: // introduce "<= b+|r|"
              this->set_r(new_index, CGAL::SMALLER);
              this->set_b(new_index, *(this->get_b()+index) + CGAL::abs(val));
              break;
            case CGAL::SMALLER:   // introduce ">=  b-|r|"
              this->set_r(new_index, CGAL::LARGER);
              this->set_b(new_index, *(this->get_b()+index) - CGAL::abs(val));
              break;
            case CGAL::EQUAL:
              if (CGAL_NTS is_positive (val)) {
                // introduce "<= b+|r|"
                this->set_r(new_index, CGAL::SMALLER);
              } else {
                // introduce ">=  b-|r|"
                this->set_r(new_index, CGAL::LARGER);
              }
              this->set_b(new_index, *(this->get_b()+index) + val);
              break;
            default:
              CGAL_qpe_assertion(false);
            }
          } else {
            this->warn1("range with identifier '%' ignored", ignored);
          }
        }

        // determine if we need to read another number:
        doItAgain = i==0 && !whitespace();
      }

      // read next token:
      t = token();
    }
    put_token_back(t);

    return true;
  }

  bool bounds_section()
  {
    std::string t = token();
    if (t != "BOUNDS") { // (Note: BOUNDS section is optional.)
      put_token_back(t);
      return true;
    }

    t = token();
    std::string bound_id;
    while (t != "QMATRIX" && t != "DMATRIX" && t != "QUADOBJ" && t != "ENDATA") {
      // process type of bound:
      enum Bound_type { LO, UP, FX, FR, MI, PL};
      Bound_type type;
      if (t=="LO")
        type = LO;
      else if (t=="UP")
        type = UP;
      else if (t=="FX")
        type = FX;
      else if (t=="FR")
        type = FR;
      else if (t=="MI")
        type = MI;
      else if (t=="PL")
        type = PL;
      else
        return
          this->err1
          ("expected 'LO', 'UP', 'FX', 'FR', 'MI', or 'PL' here but found '%'",
           t);

      // remember bound:
      const std::string bound = t;

      // find bound label; there may be several, but we only process
      // the bounds having the first bound label that occurs. This
      // label may be empty, though
      t = token(); // bound label or variable name (case of empty bound label)
      if (bound_id.size() == 0) { // first time we see a bound label / variable?
        const Index_map::const_iterator var_name = var_names.find(t);
        if (var_name != var_names.end()) { // is the token a variable?
          bound_id = " ";    // there is no bound label
          put_token_back(t); // the variable name is processed below
        } else
          bound_id = t; // we found a bound label
      } else {
        // now we already know the bound label
        if (bound_id == " ") // empty bound label?
          put_token_back(t); // the variable name is processed below
        else
          if (t != bound_id) {
            this->warn1("ignoring all bounds for bound label '%'",t);
            this->warn1("(only bounds for bound label '%' are accepted)",
                        bound_id);
          }
      }

      // find variable name;
      t = token();
      const Index_map::const_iterator var_name = var_names.find(t);
      if (var_name == var_names.end()) // unknown variable?
        return this->err1("unknown variable '%' in BOUNDS section",t);
      const unsigned int var_index = var_name->second;;

      // read value of bound, if appropriate:
      NT val;
      if (type==LO || type==UP || type==FX)
        if (!number(val))
          return this->err2("expected number after '%' in % bound",t,bound);

      // store bound:
      switch (type) {
      case FX:
        this->set_u (var_index, true, val);
        CGAL_FALLTHROUGH;
      case LO:
        this->set_l (var_index, true, val);
        break;
      case UP:
        this->set_u (var_index, true, val);
        if (val <= 0 && *(this->get_fl()+var_index)
            && *(this->get_l()+var_index) == 0)
          if (val < 0)
            this->set_l(var_index, false);
        break;
      case FR:
        this->set_u(var_index, false);
        this->set_l(var_index, false);
        break;
      case MI:
        this->set_l(var_index, false);
        break;
      case PL:
        this->set_u(var_index, false);
        break;
      default:
        CGAL_qpe_assertion(false);
      }

      // read next token:
      t = token();
    }
    put_token_back(t);

    return true;
  }

  bool qmatrix_section()
  {
    std::string t = token();
    if (t!="QMATRIX" && t!="DMATRIX" && t!="QUADOBJ") { // optional
      put_token_back(t);
      return true;
    }

    // remember section name:
    D_section = t;
    const bool multiply_by_two = t=="DMATRIX";

    t = token();
    std::string bound_id;
    while (t != "ENDATA") {
      // find first variable name;
      const Index_map::const_iterator var1_name = var_names.find(t);
      if (var1_name == var_names.end()) // unknown variable?
        return this->err2("unknown first variable '%' in '%' section",
                          t, D_section);
      const unsigned int var1_index = var1_name->second;

      // find second variable name;
      t = token();
      const Index_map::const_iterator var2_name = var_names.find(t);
      if (var2_name == var_names.end()) // unknown variable?
        return this->err2("unknown second variable '%' in '%' section",
                          t, D_section);
      const unsigned int var2_index = var2_name->second;;

      // read value:
      NT val;
      if (!number(val))
        return this->err2("expected number after '%' in section '%'",
                          t, D_section);

      // multiply by two if approriate:
      if (multiply_by_two)
        val *= NT(2);

      // set entry in D:
      int i, j;
      if (var2_index <= var1_index) {
        i = var1_index; j = var2_index;
      } else {
        i = var2_index; j = var1_index;
      }
      // rule out that we previously put a different (nonzero) value at (i,j)
      NT old_val = (*(this->get_d()+i))[j];
      if (!CGAL::is_zero(old_val) && old_val != val)
        return this->err3("nonsymmetric '%' section at variables '%' and '%'",
                          D_section, var1_name->first, var2_name->first);
      this->set_d(i, j, val);

      // read next token:
      t = token();
    }
    put_token_back(t);

    return true;
  }

};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_QP_MODELS_H
