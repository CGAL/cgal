// ============================================================================
//
// Copyright (c) 1997-2004 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/QP_solver/MPS.h
// package       : $CGAL_Package: QP_solver $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Kaspar Fischer (fischerk@inf.ethz.ch)
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: MPS Format Reader for the QP-Solver
// ============================================================================

#ifndef CGAL_QP_SOLVER_MPS_H
#define CGAL_QP_SOLVER_MPS_H

#include <istream>
#include <vector>
#include <CGAL/iterator.h>

CGAL_BEGIN_NAMESPACE

namespace QP_MPS_detail {

  template<typename V>
  struct Begin
    : public std::unary_function< V, typename V::const_iterator >
  {
    typedef typename V::const_iterator result_type;
    result_type operator () ( const V& v) const { return v.begin(); }
  };

} // QP_MPS_detail

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
	 typename ET_,  // The exact number type: this number type is
			// used to compute the rank of the matrix A,
			// if need be (see method
			// has_equalities_only_and_full_rank()).
			// Also, if you are going to invoke the
			// QP_solver using the traits
			// QP_solver_MPS_traits_d<Traits,...> then ET_
			// will be the number type the QP_solver uses
			// for its exact computations.  Requirements:
			// IT_ must be losslessly convertible to ET_.
	 typename Use_sparse_representation_for_D_=Tag_false,
                        // Use a sparse representation for D. Note:
                        // must be Tag_false currently (Tag_true not
                        // yet implemented).
	 typename Use_sparse_representation_for_A_=Tag_false>
                        // Use a sparse representation for A. Note:
                        // must be Tag_false currently (Tag_true not
                        // yet implemented).
class QP_MPS_instance {
public:
  typedef IT_ IT;
  typedef ET_ ET;
  enum Row_type { LESS_EQUAL = -1, EQUAL, GREATER_EQUAL};

private:
  typedef std::vector<IT>                 Vector;
  typedef std::vector<Vector>             Matrix;
  typedef QP_MPS_detail::Begin<Vector>    Beginner;
  typedef CGAL::Join_input_iterator_1<typename Matrix::const_iterator,
			      Beginner >  Vector_iterator;
  typedef typename Vector::const_iterator Entry_iterator;
  typedef std::vector<bool>               F_vector;
  typedef F_vector::const_iterator        F_vector_iterator;
  typedef std::vector<Row_type>           Row_type_vector;

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

  typedef typename Row_type_vector::const_iterator Row_type_iterator;
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
  unsigned int var_nr, constr_nr;
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
  Index_map var_names;

  // variables used in token() (see below):
  bool use_put_back_token;
  std::string put_back_token;

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
    
  template<typename NumberType>
  bool number(NumberType& entry) {
    whitespace();
    from >> entry;
    return from.good();
  }

  bool number(CGAL::Gmpq& entry) {
    whitespace();
    CGAL::Gmpz p,q;
    char ch;
    from >> p;
    if (!from.good()) {
      return false;
    }
    from >> ch;
    if (ch != '/') {
      return false;
    } 
    from >> q;
    if (!from.good()) {
      return false;
    }
    entry = CGAL::Gmpq(p,q);
    return true;
  }

  bool name_section();
  bool rows_section();
  bool columns_section();
  bool rhs_section();
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
  //
  // The argument strategy is a pointer to a pricing strategy (allocated
  // via operator new).  When QP_MPS_instance's destructor gets called,
  // it will call delete strategy.
  QP_MPS_instance(std::istream& in,bool use_CPLEX_convention=true,
		  int verbosity=0);

  // Returns true if and only if the instance has been properly
  // constructed (i.e., if the QP could be loaded from the MPS
  // stream).
  bool is_valid();

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
  unsigned int number_of_variables()
  {
    CGAL_qpe_assertion(is_valid());
    return var_names.size();
  }

  // Returns the number of constraints in the QP.
  //
  // Precondition: is_valid()
  unsigned int number_of_constraints()
  {
    CGAL_qpe_assertion(is_valid());
    return row_names.size();
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
  A_iterator A()
  {
    CGAL_qpe_assertion(is_valid());
    return Vector_iterator(A_.begin(),Beginner());
  }

  // Returns an iterator over the vector b (as needs to be
  // passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  B_iterator b()
  {
    CGAL_qpe_assertion(is_valid());
    return b_.begin();
  }

  // Returns an iterator over the vector c (as needs to be
  // passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  C_iterator c()
  {
    CGAL_qpe_assertion(is_valid());
    return c_.begin();
  }

  // Returns an iterator over the matrix D (as needs to be
  // passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  D_iterator D()
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
  Row_type_iterator row_types()
  {
    CGAL_qpe_assertion(is_valid());
    return row_types_.begin();
  }

  // Returns an iterator of Booleans specifying for each variable
  // whether it has a finite lower bound or not (such an iterator
  // needs to be passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  FL_iterator fl() {
    CGAL_qpe_assertion(is_valid());
    return fl_.begin();
  }

  // Returns an iterator of Booleans specifying for each variable
  // whether it has a finite upper bound or not (such an iterator
  // needs to be passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  FU_iterator fu() {
    CGAL_qpe_assertion(is_valid());
    return fu_.begin();
  }

  // Returns an iterator of IT's specifying for each variable what its
  // lower bound is (provided it is finite); such an iterator needs to
  // be passed to the constructor of class QP_solver.
  //
  // Precondition: is_valid()
  U_iterator u() {
    CGAL_qpe_assertion(is_valid());
    return u_.begin();
  }

  // Returns an iterator of IT's specifying for each variable what its
  // upper bound is (provided it is finite); such an iterator needs to
  // be passed to the constructor of class QP_solver.
  //
  // Precondition: is_valid()
  L_iterator l() {
    CGAL_qpe_assertion(is_valid());
    return l_.begin();
  }

  // Returns true iff the MPS stream could be successfully parsed and
  // the loaded QP instance has a symmetric D matrix.
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

  // Returns true iff the MPS stream could be successfully parsed and
  // the QP has only equality constraints and the coefficients (from
  // the matrix A) corresponding to these equality constraints have
  // full row rank.
  //
  // Note: This routine is expensive, more precisely O(k^3)
  // where k is the number of variables (see number_of_variables()) or
  // the number of equality constraints, whichever is larger.
  bool has_equalities_only_and_full_rank();
};

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
std::ostream& operator<<(std::ostream& o,
			 QP_MPS_instance<IT_, ET_,
			 Use_sparse_representation_for_D_,
			 Use_sparse_representation_for_A_>& qp);

template<class MPS,
	 typename Is_linear_,
	 typename Is_symmetric_,
	 typename Has_equalities_only_and_full_rank_,
	 typename Is_in_standard_form_,
	 typename IT_=typename MPS::IT,
	 typename ET_=typename MPS::ET,
	 typename D_iterator_=typename MPS::D_iterator>
class QP_solver_MPS_traits_d {
public:
  typedef IT_ IT;
  typedef ET_ ET;

public:
  typedef typename MPS::Row_type     Row_type;
  static const Row_type EQUAL =         MPS::EQUAL;
  static const Row_type LESS_EQUAL =    MPS::LESS_EQUAL;
  static const Row_type GREATER_EQUAL = MPS::GREATER_EQUAL;

  typedef typename MPS::A_iterator   A_iterator;
  typedef typename MPS::B_iterator   B_iterator;
  typedef typename MPS::C_iterator   C_iterator;
  typedef          D_iterator_       D_iterator;
  typedef typename MPS::FU_iterator  FU_iterator;
  typedef typename MPS::FL_iterator  FL_iterator;
  typedef typename MPS::U_iterator   U_iterator;
  typedef typename MPS::L_iterator   L_iterator;
  typedef typename MPS::Row_type_iterator Row_type_iterator;

  typedef Is_linear_                         Is_linear;
  typedef Is_symmetric_                      Is_symmetric;
  typedef Has_equalities_only_and_full_rank_ Has_equalities_only_and_full_rank;
  typedef Is_in_standard_form_               Is_in_standard_form;
};

CGAL_END_NAMESPACE

#include <CGAL/QP_solver/MPS.C>

#endif // CGAL_QP_SOLVER_MPS_H
