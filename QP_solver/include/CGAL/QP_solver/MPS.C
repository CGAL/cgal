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
// Author(s)     : Sven Schoenherr <sven@inf.fu-berlin.de>
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp <fransw@inf.ethz.ch>
//                 Kaspar Fischer <fischerk@inf.ethz.ch>

#include <iomanip>
#include <fstream>

#include <CGAL/QP_solver/QP_partial_filtered_pricing.h>
#include <CGAL/QP_solver/QP_partial_exact_pricing.h>

CGAL_BEGIN_NAMESPACE

// for parsing numbers
template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::number(CGAL::Gmpq& entry) {
    // accept rational or floating-point format
    std::string s = token() + " "; // routines below can't deal with EOF 
    std::istringstream from1(s), from2(s);
    return 
      number_from_quotient (entry, from1) || 
      number_from_float (entry, from2);
  }

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		     Use_sparse_representation_for_A_>::number_from_quotient(CGAL::Gmpq& entry, std::istringstream& from) {
    // reads rational in the form p/q
    //whitespace();
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

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::number_from_float(CGAL::Gmpq& entry, std::istringstream& from) {
    // reads rationals from a decimal floating-point string; 
    // e.g. "0.1" will be parsed as 1/10
    //whitespace();
    char c;

    // sign -> plus_sign
    bool plus_sign = true;
    from.get(c);
    if ((c == '+') || (c == '-'))
      plus_sign = (c == '+');
    else
      from.putback(c);
    if (!from.good()) return false;

    // digit-sequence before decimal point -> n
    CGAL::Gmpz n;
    bool digits = false; // are there digits before OR after decimal point?
    from.get(c);
    if (c != '.') {
      if (!isdigit(c)) return false;
      from.putback(c);
      from >> n;
      digits = true;
    } else {
      from.putback(c);
      n = 0;
    }
    if (!from.good()) return false;
    
    // decimal point
    bool decimal_point;
    from.get(c);
    if (c == '.') {
      decimal_point = true;
    } else {
      from.putback(c);
      decimal_point = false;
    }
    if (!from.good()) return false;  // possible decimal point is eaten

    // digit-sequence after decimal point; update n with every digit
    // found and remember according power-of-ten shift in denominator
    CGAL::Gmpz d = (plus_sign? 1 : -1);
    if (decimal_point) {
      from.get(c);
      if (!isdigit(c)) {
	from.putback(c);
      } else {
	digits = true;
	do {
	  d *= 10;
	  n = n*10 + (c-'0');
	  from.get(c);
	} while (isdigit(c));
	from.putback(c);
      }
    }
    if (!from.good()) return false; 
    
    // exponent part -> e
    int e = 0;
    from.get(c);
    if (isspace(c)) {
      from.putback(c);  // number parsed, no exponent
    } else {
      if (c == 'e' || c == 'E') {
	from >> e;      // read exponent
	if (!from.good()) return false;
      }
      from.get(c);
      if (isspace(c)) {
	from.putback(c);  // number parsed, no literal identifier
      } else {
	if ( c != 'f' && c != 'F' && c != 'l' && c != 'L') {
	  return false;    // invalid literal symbol
	} 
      }
    }
    if (!from.good()) return false;

    // now build the rational number
    if (!digits) return false; 
    // handle e
    if (e > 0) {
      while (e > 0) {
	e -= 1;
	n *= 10;
      }
    } else {
      while (e < 0) {
	e += 1;
	d *= 10;
      }
    }
    entry = CGAL::Gmpq(n,d);
    return true;
  }

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::
QP_MPS_instance(std::istream& in,bool use_CPLEX_convention,
		int verbosity)
  : verbosity_(verbosity), from(in),
    is_format_okay_(false),
    is_linear_(true),
    use_CPLEX_convention(use_CPLEX_convention),
    is_symmetric_cached(false),
    has_equalities_only_and_full_rank_cached(false),
    is_in_standard_form_cached(false),
    use_put_back_token(false)
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
  const std::string t = token();
  if (t == "RANGES") {
    err("RANGES section is not supported");
    return;
  } else
    put_token_back(t);

  // read optional BOUNDS section:
  if (!bounds_section())
    return;

  // read optional QMATRIX section:
  if (!qmatrix_section())
    return;

  // check for ENDATA:
  const std::string end = token();
  if (end != "ENDATA") {
    err1("ENDDATA expected but found '%'",end);
    return;
  }

  //std::cout << "MPS stream successfully read." << std::endl;
  is_format_okay_ = true;
}

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::is_valid()
{
  if (!is_format_okay_)
    return false;
  
  #if 0 // These tests should go into the QP-solver itself.
  // additional safety checks:
  // (Note: we do not check Has_equalities_only_and_full_rank currently,
  // as it is too expensive -- maybe add it for debugging purposes?)
  if (check_tag(Is_symmetric()) && !is_symmetric())
    return err("loaded instance does not have a symmetric D matrix but Is_symmetric is Tag_true");
  if (check_tag(Is_linear()) && !is_linear())
    return err("loaded instance is not an LP (D is nonzero) but Is_linear is Tag_true");
  if (check_tag(Is_in_standard_form()) && !is_in_standard_form())
    return err("loaded instance is not in standard form but Is_in_standard_form is Tag_true");
  #endif

  return true;
}

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
const std::string& QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::error()
{
  CGAL_qpe_assertion(!is_valid());
  return error_msg;
}

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
const std::string& QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::comment()
{
  return comment_;
}

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::is_symmetric()
{
  if (!is_format_okay_)
    return false;

  if (!is_symmetric_cached) {
    is_symmetric_ = true;
    if (!is_linear()) {
      const unsigned int var_nr = var_names.size();
      for (unsigned int i=0; i<var_nr; ++i)
	for (unsigned int j=i+1; j<var_nr; ++j)
	  if (D_[i][j] != D_[j][i]) {
	    is_symmetric_ = false; 
	    return is_symmetric_;
	  }
    }
  }
  return is_symmetric_;
}
  
template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::has_equalities_only_and_full_rank()
{
  if (!is_format_okay_)
    return false;

  if (has_equalities_only_and_full_rank_cached)
    return has_equalities_only_and_full_rank_;

  // check if we have inequalities:
  for (typename Row_type_vector::const_iterator it = row_types_.begin();
       it != row_types_.end(); ++it)
    if (*it != EQUAL) {
      has_equalities_only_and_full_rank_ = false;
      return has_equalities_only_and_full_rank_;
    }

  // perform exact (!) Gaussian Elemination to determine the rank:
  // todo: this could be made much more efficient...
  typedef Quotient<ET> QET;
  typedef std::vector<QET> V;
  typedef std::vector<V>  M;
  M A; // copy of the matrix A to work on
       // Note: A[i][j] is the element in the i-th row and j-th column.
  const int n = number_of_variables();
  const int m = number_of_constraints();
  const ET one(1);
  for (int i=0; i<m; ++i) {
    A.push_back(V());
    V &row = A.back();
    for (int j=0; j<n; ++j)
      row.push_back(QET(A_[j][i],one));
  }

  int k = 0;
  for (int j=0; j<n; ++j) {
    // Invariant: The columns 0..(j-1) of A are in echelon form and
    // column j-1 has a nonzero entry at row-index k-1.
    //
    // In this iteration of the loop we zero the entries k+1..m of
    // column j of A.  In order to do this, it might be necessary to
    // first exchange two rows.

    #if 0 // debugging code
    std::cout << "Iteration (j,k) = (" << j << "," << k << ")" << std::endl;
    for (int ii=0; ii<m; ++ii) {
      for (int jj=0; jj<n; ++jj) {
	std::cout << std::setw(15) << A[ii][jj];
      }
      std::cout << std::endl;
    }
    #endif

    // search for a suitable place k:
    bool found = false;
    int l = k;
    for (; l<m; ++l)
      if (!CGAL_NTS is_zero(A[l][j])) {
	found = true;
	break;
      }
    if (!found)
      // Here, all elements k..m of column j of A are zero, so we
      // have an echelon form (without a "step") and continue:
      continue;

    // swap rows, if necessary:
    if (k != l)
      std::swap(A[k],A[l]);
    CGAL_qpe_assertion(!CGAL_NTS is_zero(A[k][j]));

    // zero out the entries below entry k:
    for (int i=k+1; i<m; ++i) {
      // we add l times row k to row i:
      const QET l = -A[i][j]/A[k][j];
      for (int jj=0; jj<n; ++jj) // for (int jj=j+1; jj<n; ++jj)
	A[i][jj] += l*A[k][jj];
    }

    // increase echelon height:
    ++k;
  }

  has_equalities_only_and_full_rank_ = k == m;
  return has_equalities_only_and_full_rank_;
}

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::name_section()
{
  const std::string t = token();
  if (t != "NAME")
    return err("expected 'NAME'");
  name = token();
  return true;
}

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::rows_section()
{
  std::string t = token();
  if (t != "ROWS")
    return err1("expected 'ROWS' but found '%'",t);

  // read 'N', 'G', 'L', or 'E', and the name of the constraint:
  t = token();
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
	const unsigned int index = row_types_.size();
	row_types_.push_back(type == 'G'? GREATER_EQUAL :
		     (type == 'E'? EQUAL : LESS_EQUAL));
	if (row_names.find(t) != row_names.end())
	  return err1("duplicate row name '%' in section ROWS",t);
	row_names.insert(String_int_pair(t,index));
	row_by_index.push_back(t);
	b_.push_back(IT(0));
      }
      break;
    default:
      return err1(
		  "expected 'N', 'L', 'E', or 'G' in ROWS section but found '%'",
		  symbol);
    }
    t = token();
  }
  put_token_back(t);

  return true;
}

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::columns_section()
{
  std::string t = token();
  if (t != "COLUMNS")
    return err1("expected 'COLUMNS' but found '%'",t);

  t = token();
  while (t != "RHS") {
    // find variable name:
    unsigned int var_index;
    std::string col_name;
    const Index_map::const_iterator var_name = var_names.find(t);
    if (var_name == var_names.end()) { // new variable?
      var_index = var_names.size();
      col_name = t;
      var_names.insert(String_int_pair(t,var_index));
      var_by_index.push_back(t);
      A_.push_back(Vector(row_names.size(),IT(0)));
      c_.push_back(IT(0));
      fl_.push_back(true);  // default lower bound is finite...
      l_.push_back(IT(0));  // ...namely zero
      fu_.push_back(false); // default upper bound is infinite
      u_.push_back(IT());   // (dummy value) 
    } else { // variable that is already known?
      var_index = var_name->second;
      col_name = var_name->first;
    }
    //std::cout << "var is " << t << std::endl;
      
    bool doItAgain = true;
    for (int i=0; doItAgain; ++i) {
      // read row identifier:
      t = token();
      //std::cout << "row is " << t << std::endl;

      // read number:
      IT val;
      if (!number(val))
	return err2("number expected after row identifier '%' in '%' COLUMNS record",t,col_name);
      //std::cout << "val is " << val << std::endl;

      // store number:
      if (t == obj) { // objective row?
	c_[var_index] = val;
      } else { // not objective row?
	const Index_map::const_iterator row_name = row_names.find(t);
	if (row_name == row_names.end())
	  return err1("unknown row identifier '%' in section COLUMNS",t);
	A_[var_index][row_name->second] = val;
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

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::rhs_section()
{
  std::string t = token();
  if (t != "RHS")
    return err1("expected 'RHS' but found '%'",t);

  t = token();
  std::string rhs_id;
  while (t != "RANGES" && t != "BOUNDS" &&
	 t != "DMATRIX" && t != "QMATRIX" && t != "QUADOBJ" &&
	 t != "ENDATA") {
    // read rhs identifier and if it is different from the one
    // from the previous iteration, ignore the whole row:
    bool ignore;
    if (rhs_id.size() == 0) { // first time we enter the loop?
      rhs_id = t;
      ignore = false;
    } else // rhs_id already set?
      ignore =  t != rhs_id;

    bool doItAgain = true;
    for (int i=0; doItAgain; ++i) {
      // read variable identifier:
      t = token();
      //std::cout << "var is " << t << std::endl;

      // read number:
      IT val;
      if (!number(val))
	return err1("number expected after '%' in this RHS record",t);
      //std::cout << "val is " << val << std::endl;

      // store number:
      const Index_map::const_iterator row_name = row_names.find(t);
      if (row_name == row_names.end())
	return err1("unknown row identifier '%' in section RHS",t);
      if (!ignore)
	b_[row_name->second] = val;
      else {
	// todo: output warning that this particular rhs was ignored?
	//warn("
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

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::bounds_section()
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
	err1("expected 'LO', 'UP', 'FX', 'FR', 'MI', or 'PL' here but found '%'",t);
    
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
	  warn1("ignoring all bounds for bound label '%'",t);
	  warn1("(only bounds for bound label '%' are accepted)",bound_id);
	}
    }

    // find variable name;
    t = token();
    const Index_map::const_iterator var_name = var_names.find(t);
    if (var_name == var_names.end()) // unknown variable?
      return err1("unknown variable '%' in BOUNDS section",t);
    const unsigned int var_index = var_name->second;;

    // read value of bound, if appropriate:
    IT val;
    if (type==LO || type==UP || type==FX)
      if (!number(val))
	return err2("expected number after '%' in % bound",t,bound);

    // store bound:
    switch (type) {
    case FX:
      fu_[var_index] = true;
      u_ [var_index] = val;
    case LO:
      fl_[var_index] = true;
      l_ [var_index] = val;
      break;
    case UP:
      fu_[var_index] = true;
      u_ [var_index] = val;
      if (val <= 0 && fl_[var_index] == true && l_[var_index] == 0)
	if (val < 0 || !use_CPLEX_convention)
	  fl_[var_index] = false;
      break;
    case FR:
      fu_[var_index] = false;
      fl_[var_index] = false;
      break;
    case MI:
      fl_[var_index] = false;
      break;
    case PL:
      fu_[var_index] = false;
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

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_MPS_instance<IT_,ET_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::qmatrix_section()
{
  std::string t = token();
  if (t!="QMATRIX" && t!="DMATRIX" && t!="QUADOBJ") { // (Note: *MATRIX
						      // section is optional.)
    put_token_back(t);
    return true;
  }

  // remember section name:
  D_section = t;
  const bool divide_by_two = t!="DMATRIX";
  const bool only_get_lower_part = t =="QUADOBJ";

  // initialize matrix D:
  initialize_D(var_names.size(),Use_sparse_representation_for_D());

  t = token();
  std::string bound_id;
  while (t != "ENDATA") {
    // find first variable name;
    const Index_map::const_iterator var1_name = var_names.find(t);
    if (var1_name == var_names.end()) // unknown variable?
      return err1("unknown first variable '%' in D/QMATRIX section",t);
    const unsigned int var1_index = var1_name->second;;
    //std::cout << "qvar1 " << t << std::endl;
      
    // find second variable name;
    t = token();
    const Index_map::const_iterator var2_name = var_names.find(t);
    if (var2_name == var_names.end()) // unknown variable?
      return err1("unknown second variable '%' in D/QMATRIX section",t);
    const unsigned int var2_index = var2_name->second;;
    //std::cout << "qvar2 " << t << std::endl;
      
    // read value:
    IT val;
    if (!number(val))
      return err1("expected number after '%' in section QMATRIX",t);

    // divide by two if approriate:
    if (divide_by_two)
      val /= 2;

    // mark problem as nonlinear if value is nonzero:
    if (!CGAL::is_zero(val))
      is_linear_ = false;

    // set entry in D:
    set_entry_in_D(var1_index,var2_index,val,
		   Use_sparse_representation_for_D());
    if (only_get_lower_part)
      // duplicate entry if not on diagonal
      if (var1_index != var2_index) 
	set_entry_in_D(var2_index,var1_index,val,
		       Use_sparse_representation_for_D());

    // read next token:
    t = token();
  }
  put_token_back(t);

  return true;
}

template<typename IT_,
	 typename ET_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
std::ostream& operator<<(std::ostream& o,
			 QP_MPS_instance<IT_, ET_,
			 Use_sparse_representation_for_D_,
			 Use_sparse_representation_for_A_>& qp)
{
  typedef QP_MPS_instance<IT_, ET_,
    Use_sparse_representation_for_D_,
    Use_sparse_representation_for_A_> MPS;
  const unsigned int n = qp.number_of_variables();
  const unsigned int m = qp.number_of_constraints();
  
  // output general information:
  using std::endl;
  const char *yes = "yes", *no = "no";
  o << "===========" << endl
    << "MPS problem" << endl
    << "===========" << endl
    << "                       linear: "
    << (qp.is_linear()? yes : no) << endl
    << "             in standard form: "
    << (qp.is_in_standard_form()? yes : no) << endl
    << "equalities only and full rank: not checked" << endl;
  //<< (qp.has_equalities_only_and_full_rank()? yes : no) << endl;
  if (!qp.is_linear())
    o << "           symmetric D matrix: not checked" << endl
   //  << (qp.is_symmetic()? yes : no) << endl
      << "      D matrix storage format: "
      << qp.D_format_type() << endl;

  if (qp.verbosity() > 1) {
    // output c:
    o << "          number of variables: "
      << qp.number_of_variables() << endl
      << "        number of constraints: "
      << qp.number_of_constraints() << endl
      << endl
      << "objective vector: " << endl << "  ";
    std::copy(qp.c(),qp.c()+n,
	      std::ostream_iterator<IT_>(o, " "));
    o << endl;

    // output D:
    if (!qp.is_linear()) {
      o << "quadratic objective matrix: " << endl;
      typename MPS::D_iterator D = qp.D();
      for (unsigned int i=0; i<n; ++i, ++D) {
	typename MPS::D_iterator::value_type entry = *D;
	o << "  ";
	for (unsigned int j=0; j<n; ++j, ++entry)
	  o << *entry << " ";
	o << endl;
      }
    }

    // output A and b:
    o << "constraints: " << endl;
    typename MPS::B_iterator b = qp.b();
    typename MPS::A_iterator A = qp.A();
    typename MPS::Row_type_iterator r = qp.row_types();
    for (unsigned int i=0; i<m; ++i) {
      for (unsigned int j=0; j<n; ++j)
	o << "  " << A[j][i];
      if (r[i] == MPS::EQUAL)
	o << " == ";
      else if (r[i] == MPS::LESS_EQUAL)
	o << " <= ";
      else
	o << " >= ";
      o << b[i] << endl;
    }

    // output bounds:
    o << "Bounds:" << endl;
    typename MPS::FU_iterator fu = qp.fu();
    typename MPS::FL_iterator fl = qp.fl();
    typename MPS::U_iterator u = qp.u();
    typename MPS::L_iterator l = qp.l();
    for (unsigned int i=0; i<n; ++i, ++fu, ++fl, ++u, ++l) {
      o << "  x" << std::left << std::setw(3) << i << "in ";
      if (*fl)
	o << '[' << *l;
      else
	o << "(-infty";
      o << ',';
      if (*fu)
	o << *u << ']';
	else
	  o << "+infty)";
      o << endl;
    }
  }
  return o;
}

template<class MPS,
	 typename Is_linear_,
	 typename Is_symmetric_,
	 typename Has_equalities_only_and_full_rank_,
	 typename Is_in_standard_form_,
	 typename IT_,
	 typename ET_,
	 typename D_iterator_>
const typename MPS::Row_type QP_solver_MPS_traits_d<MPS,
  Is_linear_,Is_symmetric_,
  Has_equalities_only_and_full_rank_,
  Is_in_standard_form_,IT_,ET_,D_iterator_>::EQUAL;

template<class MPS,
	 typename Is_linear_,
	 typename Is_symmetric_,
	 typename Has_equalities_only_and_full_rank_,
	 typename Is_in_standard_form_,
	 typename IT_,
	 typename ET_,
	 typename D_iterator_>
const typename MPS::Row_type QP_solver_MPS_traits_d<MPS,
  Is_linear_,Is_symmetric_,
  Has_equalities_only_and_full_rank_,
  Is_in_standard_form_,IT_,ET_,D_iterator_>::LESS_EQUAL;

template<class MPS,
	 typename Is_linear_,
	 typename Is_symmetric_,
	 typename Has_equalities_only_and_full_rank_,
	 typename Is_in_standard_form_,
	 typename IT_,
	 typename ET_,
	 typename D_iterator_>
const typename MPS::Row_type QP_solver_MPS_traits_d<MPS,
  Is_linear_,Is_symmetric_,
  Has_equalities_only_and_full_rank_,
  Is_in_standard_form_,IT_,ET_,D_iterator_>::GREATER_EQUAL;

// Routines to output to MPS format:

namespace QP_MPS_detail {

  template<typename T>
  struct MPS_type_name {
    static const char *name() { return 0; }
  };
  
  template<>
  struct MPS_type_name<double> {
    static const char *name() { return "floating-point"; }
  };
  
  template<>
  struct MPS_type_name<int> {
    static const char *name() { return "integer"; }
  };
  
  template<>
  struct MPS_type_name<Gmpq> {
    static const char *name() { return "rational"; }
  };

  template<typename IT>
  struct IT_to_ET {
  };
  
  template<>
  struct IT_to_ET<double> {
    typedef Double ET;
  };
  
  template<>
  struct IT_to_ET<int> {
    typedef Gmpz ET;
  };
  
  template<>
  struct IT_to_ET<Gmpq> {
    typedef Gmpq ET;
  };

} // QP_MPS_detail

template<typename A_iterator,
	 typename B_iterator,
	 typename C_iterator,
	 typename D_iterator,
	 typename FU_iterator,
	 typename FL_iterator,
	 typename U_iterator,
	 typename L_iterator,
	 typename Row_type_iterator>
void write_MPS(std::ostream& out,
	       const std::string& number_type, // pass "" to deduce
					       // the number-type from 
                                               // U_iterator::value_type
	       const std::string& description,
	       const std::string& generator_name,
	       const std::string& problem_name,
	       const int n,
	       const int m,
	       A_iterator A,
	       B_iterator b,
	       C_iterator c,
	       D_iterator D,
	       FU_iterator fu,
	       FL_iterator fl,
	       U_iterator u,
	       L_iterator l,
	       Row_type_iterator rt)
{
  typedef typename Row_type_iterator::value_type Row_type;

  // output header:
  if (number_type.length() == 0) {
    const char *tn = QP_MPS_detail::MPS_type_name<typename U_iterator::
      value_type>::name();
    if (tn != 0)
      out << "* Number-type: " << tn << "\n";
  } else
      out << "* Number-type: " << number_type << "\n";
  out << "* Description: " << description << "\n"
      << "* Generated-by: " << generator_name << "\n"
      << "NAME " << problem_name << "\n";

  // write ROWS section:
  out << "ROWS\n"
      << "  N obj\n";                       // a row for the objective function
  for (int i=0; i<m; ++i) {
    if (rt[i] == static_cast<Row_type>(-1))
      out << "  L";
    else if (rt[i] == static_cast<Row_type>(0))
      out << "  E";
    else if (rt[i] == static_cast<Row_type>(1))
      out << "  G";
    else
      CGAL_qpe_assertion_msg(false, "incorrect row-type");
    out << " c" << i << "\n";
  }

  // output COLUMNS section:
  out << "COLUMNS\n";
  for (int i=0; i<n; ++i) {
    if (!CGAL_NTS is_zero (c[i]))
      out << "  x" << i << "  obj  " << c[i] << "\n";
    for (int j=0; j<m; ++j) { 
      if (!CGAL_NTS is_zero (A[i][j]))
	out << "  x" << i << "  c" << j << "  " << A[i][j] << "\n";
    }
  }

  // output RHS section:
  out << "RHS\n";
  for (int i=0; i<m; ++i)  
    if (!CGAL_NTS is_zero (b[i]))
      out << "  rhs c" << i << "  " << b[i] << "\n";

  // output BOUNDS section:
  out << "BOUNDS\n";
  for (int i=0; i<n; ++i) {
    #if 0 // debugging
    std::cout << "fl" << i << ": " <<  *(fl+i) << " ";
    if (*(fl+i))
      std::cout << "l=" << l[i] << " ";
    std::cout << "fu" << i << ": " <<  *(fu+i) << " ";
    if (*(fu+i))
      std::cout << "l=" << u[i] << " ";
    #endif
    if (!*(fl+i) || !CGAL::is_zero(l[i])) 
      if (*(fl+i))
	out << "  LO  BND  x" << i << "  " << l[i] << "\n";
      else
	out << "  MI  BND  x" << i << "\n";
    if (*(fu+i))
      out << "  UP  BND  x" << i << "  " << u[i] << "\n";
  }

  // check whether QMATRIX section is needed:
  bool is_linear = true;
  for (int i=0; i<n; ++i)
    for (int j=0; j<n; ++j)
      if (!CGAL::is_zero(D[i][j]))
	is_linear = false;
  
  // output QMATRIX section:
  if (!is_linear) {
    out << "QMATRIX\n";
    for (int i=0; i<n; ++i)
      for (int j=0; j<n; ++j)
	if (!CGAL::is_zero(D[i][j]))
	  out << "  x" << i << "  x" << j << "  " << 2*D[i][j] << "\n";
  }    

  // output end:
  out << "ENDATA\n";
}

CGAL_END_NAMESPACE
