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
// $URL: svn+ssh://gaertner@scm.gforge.inria.fr/svn/cgal/trunk/QP_solver/include/CGAL/QP_models.h $
// $Id: QP_models.h 33922 2006-09-05 12:32:25Z gaertner $
// 
//
// Author(s)     : Bernd Gaertner <gaertner@inf.ethz.ch>

#ifndef CGAL_QP_MODELS_H
#define CGAL_QP_MODELS_H

#include <CGAL/basic.h>
#include <CGAL/iterator.h>
#include <vector> 
#include <iomanip>
#include <istream>
#include <sstream>

#define QP_MODEL_ITERATOR_TYPES \
  typedef typename Base::A_iterator A_iterator;\
  typedef typename Base::B_iterator B_iterator;\
  typedef typename Base::R_iterator R_iterator;\
  typedef typename Base::FL_iterator FL_iterator;\
  typedef typename Base::L_iterator L_iterator;\
  typedef typename Base::FU_iterator FU_iterator;\
  typedef typename Base::U_iterator U_iterator;\
  typedef typename Base::D_iterator D_iterator; \
  typedef typename Base::C_iterator C_iterator; \
  typedef typename Base::value_type value_type;
// end QP_MODEL_ITERATOR_TYPES


CGAL_BEGIN_NAMESPACE

//  default iterator types to be used in LP / Nonnegative models
template <class Iterator>
class QP_Default {
private:
  typedef typename std::iterator_traits<Iterator>::value_type value_type;
public:
  typedef Const_oneset_iterator<value_type> 
  It_1d; // 1-dimensional random access iterator for a constant value
  typedef Const_oneset_iterator<Const_oneset_iterator<value_type> > 
  It_2d; // 2-dimensional random access iterator for a constant value
};


// models of QuadraticProgram
// ==========================

// QP_from_iterators
// -----------------
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
class QP_from_iterators 
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
  typedef typename std::iterator_traits<C_iterator>::value_type value_type;
private:
  // data
  const int n_;
  const int m_;
  const A_iterator a_it;
  const B_iterator b_it; 
  const R_iterator r_it;
  const FL_iterator fl_it;
  const L_iterator l_it;
  const FU_iterator fu_it;
  const U_iterator u_it; 
  const D_iterator d_it;
  const C_iterator c_it;
  const value_type c_0; // constant term in objective function
 
public:
  // construction
  QP_from_iterators (
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
		     const value_type& c0 = value_type(0)
		     )
    : n_ (n), m_ (m), a_it (a), b_it (b), r_it (r), fl_it (fl), l_it (l), 
      fu_it (fu), u_it (u), d_it (d), c_it (c), c_0 (c0)    
  {}

  // access
  int n() const {return n_;}
  int m() const {return m_;}
  const A_iterator& a() const {return a_it;}
  const B_iterator& b() const {return b_it;}
  const R_iterator& r() const {return r_it;}  
  const FL_iterator& fl() const {return fl_it;}
  const L_iterator& l() const {return l_it;}
  const FU_iterator& fu() const {return fu_it;}
  const U_iterator& u() const {return u_it;}
  const D_iterator& d() const {return d_it;}
  const C_iterator& c() const {return c_it;}
  const value_type& c0() const {return c_0;}
};

// global function make_QP_from_iterators
// --------------------------------------
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
QP_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, D_it, C_it>
make_QP_from_iterators (
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
   typename std::iterator_traits<C_it>::value_type c0)
{
  return QP_from_iterators
    <A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, D_it, C_it>
    (n, m, a, b, r, fl, l, fu, u, d, c, c0);
}	

// QP_from_pointers
// ----------------
template <typename NT>
class QP_from_pointers : 
  public QP_from_iterators <NT**, NT*, CGAL::Comparison_result*, 
			    bool*, NT*, bool*, NT*, NT**, NT*>
{
private:
  typedef QP_from_iterators <NT**, NT*, CGAL::Comparison_result*, 
			     bool*, NT*, bool*, NT*, NT**, NT*> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  QP_from_pointers (
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
		     const value_type& c0 = value_type(0)
		     )
    : Base (n, m, a, b, r, fl, l, fu, u, d, c, c0)
  {}  
};

// LP_from_iterators
// -----------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename FL_it,  // for finiteness of lower bounds (value type bool)
  typename L_it,   // for lower bounds
  typename FU_it,  // for finiteness of upper bounds (value type bool)
  typename U_it,   // for upper bounds
  typename C_it >  // for objective function c
class LP_from_iterators : 
  public QP_from_iterators <A_it, B_it, R_it, FL_it, L_it, FU_it, U_it,
			    typename QP_Default<A_it>::It_2d, C_it>
{
private:
  typedef QP_from_iterators <A_it, B_it, R_it, FL_it, L_it, FU_it, U_it,
			     typename QP_Default<B_it>::It_2d, C_it> Base;
  typedef typename QP_Default<B_it>::It_2d Const_D_iterator;
public:
   QP_MODEL_ITERATOR_TYPES;
   LP_from_iterators (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,
		     const FL_iterator& fl,
		     const L_iterator& l,
		     const FU_iterator& fu,
		     const U_iterator& u,
		     const C_iterator& c,
		     const value_type& c0 = value_type(0)
		     )
    : Base (n, m, a, b, r, fl, l, fu, u, 
	    Const_D_iterator(value_type(0)), c, c0)
  {}  
}; 

// LP_from_pointers
// ----------------
template <typename NT>
class LP_from_pointers : 
  public LP_from_iterators <NT**, NT*, CGAL::Comparison_result*, bool*, 
			    NT*, bool*, NT*, NT*>
{
private:
  typedef  LP_from_iterators <NT**, NT*, CGAL::Comparison_result*, bool*, 
			    NT*, bool*, NT*, NT*> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  LP_from_pointers (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,
		     const FL_iterator& fl,
		     const L_iterator& l,
		     const FU_iterator& fu,
		     const U_iterator& u,
		     const C_iterator& c,
		     const value_type& c0 = value_type(0)
		     )
    : Base (n, m, a, b, r, fl, l, fu, u, c, c0)
  {}  
};

// Nonnegative_QP_from_iterators
// -----------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename D_it,   // for quadratic matrix D (rowwise)
  typename C_it >  // for objective function c
class Nonnegative_QP_from_iterators : 
  public QP_from_iterators <A_it, B_it, R_it, 
			    typename QP_Default<bool*>::It_1d, 
			    typename QP_Default<B_it>::It_1d,
			    typename QP_Default<bool*>::It_1d, 
			    typename QP_Default<B_it>::It_1d,
			    D_it, C_it>
{
private:
  typedef  QP_from_iterators <A_it, B_it, R_it, 
			      typename QP_Default<bool*>::It_1d, 
			      typename QP_Default<B_it>::It_1d,
			      typename QP_Default<bool*>::It_1d, 
			      typename QP_Default<B_it>::It_1d,
			      D_it, C_it> Base;
  typedef typename QP_Default<bool*>::It_1d Const_FLU_iterator;
  typedef typename QP_Default<B_it>::It_1d Const_LU_iterator;
public:
   QP_MODEL_ITERATOR_TYPES;
   Nonnegative_QP_from_iterators (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,		     
		     const D_iterator& d,
		     const C_iterator& c,
		     const value_type& c0 = value_type(0)
		     )
    : Base (n, m, a, b, r, 
	    Const_FLU_iterator(true), Const_LU_iterator(value_type(0)), 
	    Const_FLU_iterator(false), Const_LU_iterator(value_type(0)), 
	    d, c, c0)
  {}  
};

// Nonnegative_QP_from_pointers
// ----------------------------
template <typename NT>
class Nonnegative_QP_from_pointers : 
  public Nonnegative_QP_from_iterators <NT**, NT*, CGAL::Comparison_result*, 
					NT**, NT*>
{
private:
  typedef Nonnegative_QP_from_iterators <NT**, NT*, CGAL::Comparison_result*, 
					NT**, NT*> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  Nonnegative_QP_from_pointers (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,
		     const D_iterator& d,
		     const C_iterator& c,
		     const value_type& c0 = value_type(0)
		     )
    : Base (n, m, a, b, r, d, c, c0)
  {}  
};


 
// Nonnegative_LP_from_iterators
// -----------------------------
template <
  typename A_it,   // for constraint matrix A (columnwise)
  typename B_it,   // for right-hand side b 
  typename R_it,   // for relations (value type Comparison)
  typename C_it >  // for objective function c
class Nonnegative_LP_from_iterators : 
  public QP_from_iterators <A_it, B_it, R_it, 
			    typename QP_Default<bool*>::It_1d, 
			    typename QP_Default<B_it>::It_1d,
			    typename QP_Default<bool*>::It_1d, 
			    typename QP_Default<B_it>::It_1d,
			    typename QP_Default<B_it>::It_2d, C_it>
{
private:
  typedef  QP_from_iterators <A_it, B_it, R_it, 
			      typename QP_Default<bool*>::It_1d, 
			      typename QP_Default<B_it>::It_1d,
			      typename QP_Default<bool*>::It_1d, 
			      typename QP_Default<B_it>::It_1d,
			      typename QP_Default<B_it>::It_2d, C_it> Base;
  typedef typename QP_Default<bool*>::It_1d Const_FLU_iterator;
  typedef typename QP_Default<B_it>::It_1d Const_LU_iterator;
  typedef typename QP_Default<B_it>::It_2d Const_D_iterator;
public:
   QP_MODEL_ITERATOR_TYPES;
   Nonnegative_LP_from_iterators (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,		     
		     const C_iterator& c,
		     const value_type& c0 = value_type(0)
		     )
    : Base (n, m, a, b, r, 
	    Const_FLU_iterator(true), Const_LU_iterator(value_type(0)), 
	    Const_FLU_iterator(false), Const_LU_iterator(value_type(0)), 
	    Const_D_iterator(value_type(0)), c, c0)
  {}  
}; 

// Nonnegative_LP_from_pointers
// ----------------------------
template <typename NT>
class Nonnegative_LP_from_pointers : 
  public Nonnegative_LP_from_iterators <NT**, NT*, CGAL::Comparison_result*, 
					NT*>
{
private:
  typedef Nonnegative_LP_from_iterators <NT**, NT*, CGAL::Comparison_result*, 
					NT*> Base;
public:
  QP_MODEL_ITERATOR_TYPES;
  Nonnegative_LP_from_pointers (
		     int n, int m, // number of variables / constraints
		     const A_iterator& a, 
		     const B_iterator& b,
		     const R_iterator& r,
		     const C_iterator& c,
		     const value_type& c0 = value_type(0)
		     )
    : Base (a, b, r, c, c0)
  {}  
};

// QP_from_mps
// -----------
namespace QP_from_mps_detail {

  template<typename V>
  struct Begin
    : public std::unary_function< V, typename V::const_iterator >
  {
    typedef typename V::const_iterator result_type;
    result_type operator () ( const V& v) const { return v.begin(); }
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
	 typename Use_sparse_representation_for_D_=Tag_false,
                        // Use a sparse representation for D. Note:
                        // must be Tag_false currently (Tag_true not
                        // yet implemented).
	 typename Use_sparse_representation_for_A_=Tag_false>
                        // Use a sparse representation for A. Note:
                        // must be Tag_false currently (Tag_true not
                        // yet implemented).
class QP_from_mps {
public:
  typedef IT_ IT;

public: // undocumented types, should be considered private:
  typedef std::vector<IT>                 Vector;
  typedef std::vector<Vector>             Matrix;
  typedef QP_from_mps_detail::Begin<Vector>    Beginner;
  typedef CGAL::Join_input_iterator_1<typename Matrix::const_iterator,
			      Beginner >  Vector_iterator;
  typedef typename Vector::const_iterator Entry_iterator;
  typedef std::vector<bool>               F_vector;
  typedef F_vector::const_iterator        F_vector_iterator;
  typedef std::vector<CGAL::Comparison_result> Row_type_vector;

public:
  // iterators over the input matrices and vectors: 
  typedef Vector_iterator   A_iterator;
  typedef Entry_iterator    B_iterator;
  typedef Entry_iterator    C_iterator;
  typedef Vector_iterator   D_iterator;
  typedef Const_oneset_iterator< Const_oneset_iterator<IT> >
                            Zero_D_iterator;
  typedef F_vector_iterator FU_iterator;
  typedef F_vector_iterator FL_iterator;
  typedef Entry_iterator    U_iterator;
  typedef Entry_iterator    L_iterator;
  typedef IT                value_type;

  typedef typename Row_type_vector::const_iterator R_iterator;
  typedef Use_sparse_representation_for_D_   Use_sparse_representation_for_D;
  typedef Use_sparse_representation_for_A_   Use_sparse_representation_for_A;

private:
  typedef std::pair<std::string,unsigned int> String_int_pair;
  typedef std::map<std::string,unsigned int> Index_map;

private:
  const int verbosity_;
  std::istream& from;
  std::string error_msg;
  bool is_format_okay_;
  bool is_linear_;
  const bool use_CPLEX_convention;
  IT c_0;                      // constant term in objective function
  Vector b_, c_;              // vectors b and c
  Matrix A_, D_;              // matrices A and D
  std::string D_section;      // name of the section from which D was read
  Row_type_vector row_types_; // equality type for each row
  Vector u_, l_;              // upper and lower bounds
  F_vector fu_, fl_;          // whether the lower/upper bound is finite or not

  // cached data:
  bool is_symmetric_cached, is_symmetric_;
  bool has_equalities_only_and_full_rank_cached,
    has_equalities_only_and_full_rank_;
  bool is_in_standard_form_cached, is_in_standard_form_;

  // data from MPS file:
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

public: // unofficial helper routines used in master_mps_to_derivatives.C;
        // these should be considered private for CGAL end-users:

  const Matrix& A_matrix() { return A_; }
  const Vector& b_vector() { return b_; }
  const Row_type_vector& row_types_vector() { return row_types_; }

private: // helpers to deal with sparse/dense representation:

  void initialize_D(unsigned int dimension,const Tag_true) // sparse case
  {
    CGAL_qpe_assertion_msg(false, "not implemented yet");
  }

  void initialize_D(unsigned int dimension,const Tag_false) // dense case
  {
    for (unsigned int i=0; i<dimension; ++i)
      D_.push_back(Vector(dimension,IT(0)));
  }

  void set_entry_in_D(unsigned int i,unsigned int j,const IT& val,
		      const Tag_true)                       // sparse case
  {
    CGAL_qpe_assertion_msg(false, "not implemented yet");
  }

  void set_entry_in_D(unsigned int i,unsigned int j,const IT& val,
		      const Tag_false)                      // dense case
  {
    D_[i][j] = val;
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
    
  // here is a general template for halving a number; we must
  // have a specialization for Gmpz that performs exact division
  // (and therefore fails if the result is incorrect)
  template<typename NumberType>
  void halve(NumberType& entry) {
    entry /= 2;
  }

  void halve(CGAL::Gmpz& entry) {
    entry = CGAL::exact_division(entry,2);    
  }

  // here is the general template, and the C-file contains
  // implementations of a specialization for Gmpq
  template<typename NumberType>
  bool number(NumberType& entry) {
    whitespace();
    from >> entry;
    return from.good();
  }
  
  bool number(CGAL::Gmpq& entry) {
    // accept rational or floating-point format
    std::string s = token() + " "; // routines below can't deal with EOF 
    std::istringstream from1(s), from2(s);
    return 
      number_from_quotient (entry, from1) || 
      number_from_float (entry, from2);
  }
  
  bool number_from_quotient(CGAL::Gmpq& entry, std::istringstream& from);
  bool number_from_float(CGAL::Gmpq& entry, std::istringstream& from);

  bool name_section();
  bool rows_section();
  bool columns_section();
  bool rhs_section();
  bool ranges_section();
  bool bounds_section();
  bool qmatrix_section();

private:
  D_iterator D(const Tag_true);
  D_iterator D(const Tag_false);

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
  unsigned int n() const
  {
    CGAL_qpe_assertion(is_valid());
    return var_names.size();
  }

  // Returns the number of constraints in the QP.
  //
  // Precondition: is_valid()
  unsigned int m() const
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
  const std::string& name_of_variable(unsigned int i)
  {
    CGAL_qpe_assertion(0<=i && i<n());
    return var_by_index[i];
//     for (Index_map::const_iterator it = var_names.begin();
// 	 it != var_names.end();
// 	 ++it)
//       if ((*it).second == i)
//         return (*it).first;

//     CGAL_qpe_assertion(false); // should never get here; following
// 			       // only to fix compiler warnings
//     static std::string dummy;
//     return dummy;
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

  // Initializes the internal D matrix to the zero matrix so that D()
  // will return an iterator over a zero matrix.  Normally, there
  // should not be any need to call this routine; it is used by the
  // QP_solver testsuite to solve an LP as a QP, and for this we need
  // a zero D matrix.
  //
  // Note: If Use_sparse_representation_for_D_ is Tag_false, this
  // routine allocates a zero (n x n)-matrix (where n is
  // number_of_variables()).  If you known in advance that you need a
  // zero D matrix, you might want to use zero_D() (see below).
  //
  // Precondition: is_linear()
  void make_zero_D() 
  {
    CGAL_qpe_assertion(is_linear());
    initialize_D(var_names.size(),Use_sparse_representation_for_D());
  }

  // Returns an iterator over the matrix A (as needs to be
  // passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  A_iterator a() const
  {
    CGAL_qpe_assertion(is_valid());
    return Vector_iterator(A_.begin(),Beginner());
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
  D_iterator d() const
  {
    CGAL_qpe_assertion(is_valid());
    return Vector_iterator(D_.begin(),Beginner());
  }

  // Returns an iterator over a matrix D that is zero (as needs to be
  // passed to the constructor of class QP_solver if you for instance
  // want to solve an LP as a QP with a zero D-matrix for test
  // purposes).
  //
  // Precondition: is_valid()
  Zero_D_iterator zero_D()
  {
    CGAL_qpe_assertion(is_valid());
    return Zero_D_iterator(Const_oneset_iterator<IT>(0));
  }

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
  // the loaded QP instance has a symmetric D matrix (if it has a D
  // matrix at all).
  bool is_symmetric();

  // Returns true iff the MPS stream could be successfully parsed and
  // the loaded QP instance has a D matrix that is zero (i.e., iff the
  // QP is actually an LP).
  bool is_linear()
  {
    return is_format_okay_ && is_linear_;
  }

  // Returns true iff the MPS stream could be successfully parsed and
  // the loaded QP instance is in standard form.
  bool is_in_standard_form()
  {
    if (!is_format_okay_)
      return false;

    if (!is_in_standard_form_cached) {
      for (unsigned int i=0; i<var_names.size(); ++i)
	if (fl_[i] == false || l_[i] != 0 ||
	    fu_[i] == true) {
	  is_in_standard_form_ = false;
	  return is_in_standard_form_;
	}
      is_in_standard_form_ = true;
    }
    return is_in_standard_form_;
  }
};

template<typename IT_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
std::ostream& operator<<(std::ostream& o,
			 QP_from_mps<IT_,
			 Use_sparse_representation_for_D_,
			 Use_sparse_representation_for_A_>& qp);


CGAL_END_NAMESPACE

#include <CGAL/QP_solver/QP_models_impl.h>

#endif // CGAL_QP_MODELS_H
