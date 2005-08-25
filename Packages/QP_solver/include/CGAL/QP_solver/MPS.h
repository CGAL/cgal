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

template<class Traits>
class QP_MPS_instance;

template<typename IT_,  // Warning: IT_ must support EXACT multiplication
                        // by 2 (so x*2 must be exactly representable in IT_).
	                // The reason is that the MPS format stores the matrix
	                // D/2 and the QP-solver needs D, so we need to
                        // multiply by 2.
	 typename ET_,
	 typename Is_linear_,
	 typename Is_symmetric_,
	 typename Has_equalities_only_and_full_rank_,
	 typename Is_in_standard_form_,
	 typename Use_sparse_representation_for_D_=Tag_false> // todo: maybe
                                                    // add optional sparse
                                                    // representation for A
                                                    // as well?
                                                    // Note: must be Tag_false
                                                    // currently (Tag_true
                                                    // not yet implemented).
class QP_solver_MPS_traits_d {
 template<class Traits> friend class QP_MPS_instance;

public:
  typedef IT_ IT;                       // input type (i.e., type of the
                                        // numbers in the MPS stream)
  typedef ET_ ET;                       // exact number type

private:
  typedef std::vector<IT>     Vector;
  typedef std::vector<Vector> Matrix;
  typedef QP_MPS_detail::Begin<Vector> Beginner;
  typedef CGAL::Join_input_iterator_1<typename Matrix::const_iterator,
			      Beginner >
                                          Vector_iterator;
  typedef typename Vector::const_iterator Entry_iterator;
  typedef std::vector<bool> F_vector;
  typedef F_vector::const_iterator F_vector_iterator;

public:
  enum Row_type { LESS_EQUAL = -1, EQUAL, GREATER_EQUAL};

private:
  typedef std::vector<Row_type> Row_type_vector;

public:
  // iterators over the input matrices and vectors:
  typedef Vector_iterator   A_iterator;
  typedef Entry_iterator    B_iterator;
  typedef Entry_iterator    C_iterator;
  typedef Vector_iterator   D_iterator;
  typedef F_vector_iterator FU_iterator;
  typedef F_vector_iterator FL_iterator;
  typedef Vector_iterator   U_iterator;
  typedef Vector_iterator   L_iterator;

  typedef typename Row_type_vector::const_iterator Row_type_iterator;
  typedef Is_linear_                         Is_linear;
  typedef Is_symmetric_                      Is_symmetric;
  typedef Has_equalities_only_and_full_rank_ Has_equalities_only_and_full_rank;
  typedef Is_in_standard_form_               Is_in_standard_form;

  typedef Use_sparse_representation_for_D_   Use_sparse_representation_for_D;
};

template<class Traits>
class QP_MPS_instance {
private:
  typedef typename Traits::IT IT;
  typedef typename Traits::ET ET;
  typedef typename Traits::Vector Vector;
  typedef typename Traits::F_vector F_vector;
  typedef typename Traits::Matrix Matrix;
  typedef typename Traits::Entry_iterator Entry_iterator;
  typedef typename Traits::Vector_iterator Vector_iterator;
  typedef typename Traits::Beginner Beginner;
  typedef typename Traits::Row_type Row_type;
  typedef typename Traits::Row_type_vector Row_type_vector;
  typedef typename Traits::Row_type_iterator Row_type_iterator;
  typedef typename Traits::Is_linear Is_linear;
  typedef typename Traits::Is_symmetric Is_symmetric;
  typedef typename Traits::Has_equalities_only_and_full_rank 
                           Has_equalities_only_and_full_rank;
  typedef typename Traits::Is_in_standard_form Is_in_standard_form;

  typedef typename Traits::Use_sparse_representation_for_D
                           Use_sparse_representation_for_D;

  typedef typename Traits::A_iterator  A_iterator;
  typedef typename Traits::B_iterator  B_iterator;
  typedef typename Traits::C_iterator  C_iterator;
  typedef typename Traits::D_iterator  D_iterator;
  typedef typename Traits::FU_iterator FU_iterator;
  typedef typename Traits::FL_iterator FL_iterator;
  typedef typename Traits::U_iterator  U_iterator;
  typedef typename Traits::L_iterator  L_iterator;

  typedef std::pair<std::string,unsigned int> String_int_pair;
  typedef std::map<std::string,unsigned int> Index_map;

private:
  std::istream& from;
  std::string error_msg;
  bool is_format_okay_;
  bool is_linear_;
  const bool use_CPLEX_convention;
  unsigned int var_nr, constr_nr;
  Vector b_, c_;              // vectors b and c
  Matrix A_, D_;              // matrices A and D
  Row_type_vector row_types_; // equality type for each row
  Vector u_, l_;              // upper and lower bounds
  F_vector fu_, fl_;          // whether the lower/upper bound is finite or not
  QP_pricing_strategy<Traits> *strategy;

  // data from MPS file:
  std::string name;  // from the NAME section
  std::string obj;   // name of the objective "constraint"
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
		      const Tag_true) // sparse case
  {
    CGAL_qpe_assertion_msg(false, "not implemented yet");
  }

  void set_entry_in_D(unsigned int i,unsigned int j,const IT& val,
		      const Tag_false) // dense case
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
    if (use_put_back_token)
      return false;
    bool lineBreakFound = false;
    char c;
    while (from.get(c)) {
      if (!isspace(c)) {
	if (c!='$' && c!='*')
	  from.putback(c);
	else { // comment?
	  while (from.get(c) && c!='\r' && c!='\n')
	    ; // nop
	  lineBreakFound = true;
	}
	break;
      }
      if (c=='\r' || c=='\n')
	lineBreakFound = true;
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
    from >> token;
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
		  QP_pricing_strategy<Traits> *strategy = 0);

  // Destructor.
  ~QP_MPS_instance();

  // Returns true if and only if the instance has been properly
  // constructed (i.e., if the QP could be loaded from the MPS
  // stream).
  bool is_valid();

  // If is_valid() returns false, this routine returns an error
  // string describing the problem that was encountered.
  //
  // Precondition: !is_valid()
  const std::string error();
  
  // Returns the number of variables in the QP.
  //
  // Precondition: is_valid()
  unsigned int number_of_variables();

  // Returns the number of constraints in the QP.
  //
  // Precondition: is_valid()
  unsigned int number_of_constraints();

  // Returns an iterator over the matrix A (as needs to be
  // passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  A_iterator A();

  // Returns an iterator over the vector b (as needs to be
  // passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  B_iterator b();

  // Returns an iterator over the vector c (as needs to be
  // passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  C_iterator c();

  // Returns an iterator over the matrix D (as needs to be
  // passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  D_iterator D();

  // Returns an iterator (of value-type Row_type) over the constraint
  // types (as needs to be passed to the constructor of class
  // QP_solver).
  //
  // Precondition: is_valid()
  Row_type_iterator row_types();

  // Returns an iterator of Booleans specifying for each variable
  // whether it has a finite lower bound or not (such an iterator
  // needs to be passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  FL_iterator fl() {
    CGAL_qpe_assertion(!is_valid());
    return fl_.begin();
  }

  // Returns an iterator of Booleans specifying for each variable
  // whether it has a finite upper bound or not (such an iterator
  // needs to be passed to the constructor of class QP_solver).
  //
  // Precondition: is_valid()
  FU_iterator fu() {
    CGAL_qpe_assertion(!is_valid());
    return fu_.begin();
  }

  // Returns an iterator of IT's specifying for each variable what its
  // lower bound is (provided it is finite); such an iterator needs to
  // be passed to the constructor of class QP_solver.
  //
  // Precondition: is_valid()
  FU_iterator u() {
    CGAL_qpe_assertion(!is_valid());
    return u_.begin();
  }

  // Returns an iterator of IT's specifying for each variable what its
  // upper bound is (provided it is finite); such an iterator needs to
  // be passed to the constructor of class QP_solver.
  //
  // Precondition: is_valid()
  FU_iterator l() {
    CGAL_qpe_assertion(!is_valid());
    return l_.begin();
  }

  // Returns the default pricing strategy for solving this QP.
  //
  // Precondition: is_valid()
  QP_pricing_strategy<Traits>& default_pricing_strategy() {
    if (strategy == 0)
      strategy = new QP_partial_filtered_pricing<Traits>;
    return *strategy;
  }
  

  // Returns true iff the loaded QP instance has a symmetric D matrix.
  //
  // Precondition: is_valid()
  bool is_symmetric();

  // Returns true iff the loaded QP instance has a D matrix that is zero
  // (i.e., iff the QP is actually an LP).
  //
  // Precondition: is_valid()
  bool is_linear()
  {
    return is_linear_;
  }

  // Returns true iff the loaded QP instance is in standard form.
  //
  // Precondition: is_valid()
  bool is_in_standard_form()
  {
    for (unsigned int i=0; i<var_names.size(); ++i)
      if (fl_[i] == false || l_[i] != 0 ||
	  fu_[i] == true)
	return false;
    return true;
  }
};

template<class Traits>
std::ostream& operator<<(std::ostream& o,QP_MPS_instance<Traits>& qp);

CGAL_END_NAMESPACE

#include <CGAL/QP_solver/MPS.C>

#endif // CGAL_QP_SOLVER_MPS_H
