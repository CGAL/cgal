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
// $URL: svn+ssh://gaertner@scm.gforge.inria.fr/svn/cgal/trunk/QP_solver/include/CGAL/QP_solver/QP_models_impl.h $
// $Id: QP_models_impl.h 33922 2006-09-05 12:32:25Z gaertner $
// 
//
// Author(s)     : Sven Schoenherr <sven@inf.fu-berlin.de>
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp <fransw@inf.ethz.ch>
//                 Kaspar Fischer <fischerk@inf.ethz.ch>

#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/QP_functions.h>

CGAL_BEGIN_NAMESPACE

// for parsing rational numbers
template<typename IT_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_from_mps<IT_,
		Use_sparse_representation_for_D_,
		     Use_sparse_representation_for_A_>::number_from_quotient(CGAL::Gmpq& entry, std::istringstream& from) {
    // reads rational in the form p/q
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
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_from_mps<IT_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::number_from_float(CGAL::Gmpq& entry, std::istringstream& from) {
    // reads rationals from a decimal floating-point string; 
    // e.g. "0.1" will be parsed as 1/10
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
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
QP_from_mps<IT_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::
QP_from_mps(std::istream& in,bool use_CPLEX_convention,
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
    err1("ENDDATA expected but found '%'",end);
    return;
  }

  //std::cout << "MPS stream successfully read." << std::endl;
  is_format_okay_ = true;
}

template<typename IT_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_from_mps<IT_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::is_valid() const
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
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
const std::string& QP_from_mps<IT_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::error()
{
  CGAL_qpe_assertion(!is_valid());
  return error_msg;
}

template<typename IT_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
const std::string& QP_from_mps<IT_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::comment()
{
  return comment_;
}

template<typename IT_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_from_mps<IT_,
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
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_from_mps<IT_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::name_section()
{
  const std::string t = token();
  if (t != "NAME")
    return err("expected 'NAME'");
  // NAME: everything found until line break; whitespaces are allowed
  char c;
  std::string token;
  std::string whitespaces;
  if (whitespace()) 
    // line break eaten, name is empty
    return true;
  do {
    from.get(c);
    if (c == '\r' || c == '\n') return true;
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

template<typename IT_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_from_mps<IT_,
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
	row_types_.push_back(type == 'G'? CGAL::LARGER :
		     (type == 'E'? CGAL::EQUAL : CGAL::SMALLER));
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
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_from_mps<IT_,
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
      
    bool doItAgain = true;
    for (int i=0; doItAgain; ++i) {
      // read row identifier:
      t = token();

      // read number:
      IT val;
      if (!number(val))
	return err2("number expected after row identifier '%' in '%' COLUMNS record",t,col_name);

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
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_from_mps<IT_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::rhs_section()
{
  c_0 = IT(0);  // no constant term yet
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
      IT val;
      if (!number(val))
	return err1("number expected after '%' in this RHS record",t);

      // store number:
      const Index_map::const_iterator row_name = row_names.find(t);
      if (row_name == row_names.end()) {
	// no corresponding constraint; is it the constant term?
	if (t == obj) 
	  c_0 = -val;
	else 
	  return err1("unknown row identifier '%' in section RHS",t);
      } else {
	// we have an actual constraint
	if (!ignore) {
	  b_[row_name->second] = val;
	} else {
	  warn1("rhs with identifier '%' ignored", ignored);
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

template<typename IT_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_from_mps<IT_,
		Use_sparse_representation_for_D_,
		Use_sparse_representation_for_A_>::ranges_section()
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
      IT val;
      if (!number(val))
	return err1("number expected after '%' in this RANGES record",t);

      // duplicate the constraint, depending on sign of val and type
      // of constraint
      const Index_map::const_iterator row_name = row_names.find(t);
      if (row_name == row_names.end()) {
	  return err1("unknown row identifier '%' in section RANGES",t);
      } else {
	if (!ignore) {
	  int index = row_name->second;
	  CGAL::Comparison_result type = row_types_[index];
	  // duplicate the row, unless it has already been duplicated
	  const Index_map::const_iterator duplicated_row_name = 
	    duplicated_row_names.find(t);
	  if (duplicated_row_name != duplicated_row_names.end())
	    return err1("duplicate row identifier '%' in section RANGES",t);
	  duplicated_row_names.insert(*row_name);
	  for (unsigned int j=0; j<var_names.size(); ++j) {
	    A_[j].push_back(A_[j][index]);
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
	    row_types_.push_back(CGAL::SMALLER);
	    b_.push_back(b_[index] + CGAL::abs(val));
	    break;
	  case CGAL::SMALLER:   // introduce ">=  b-|r|" 
	    row_types_.push_back(CGAL::LARGER);
	    b_.push_back(b_[index] - CGAL::abs(val));
	    break;
	  case CGAL::EQUAL:   
	    if (CGAL_NTS is_positive (val)) {
	      // introduce "<= b+|r|"
	      row_types_.push_back(CGAL::SMALLER);
	    } else {
	      // introduce ">=  b-|r|"  
	      row_types_.push_back(CGAL::LARGER);  
	    }
	    b_.push_back(b_[index] + val);
	    break;
	  default:
	    CGAL_qpe_assertion(false);
	  }  
	} else {
	  warn1("range with identifier '%' ignored", ignored);
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
	 


template<typename IT_,
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_from_mps<IT_,
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
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
bool QP_from_mps<IT_,
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
      halve (val);

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
	 typename Use_sparse_representation_for_D_,
	 typename Use_sparse_representation_for_A_>
std::ostream& operator<<(std::ostream& o,
			 QP_from_mps<IT_, 
			 Use_sparse_representation_for_D_,
			 Use_sparse_representation_for_A_>& qp)
{
  typedef QP_from_mps<IT_, 
    Use_sparse_representation_for_D_,
    Use_sparse_representation_for_A_> MPS;
  const unsigned int n = qp.n();
  const unsigned int m = qp.m();
  
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
      << qp.n() << endl
      << "        number of constraints: "
      << qp.m() << endl
      << endl
      << "objective vector: " << endl << "  ";
    std::copy(qp.c(),qp.c()+n,
	      std::ostream_iterator<IT_>(o, " "));
    o << endl;

    // output D:
    if (!qp.is_linear()) {
      o << "quadratic objective matrix: " << endl;
      typename MPS::D_iterator D = qp.d();
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
    typename MPS::A_iterator A = qp.a();
    typename MPS::R_iterator r = qp.r();
    for (unsigned int i=0; i<m; ++i) {
      for (unsigned int j=0; j<n; ++j)
	o << "  " << A[j][i];
      if (r[i] == CGAL::EQUAL)
	o << " == ";
      else if (r[i] == CGAL::SMALLER)
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

CGAL_END_NAMESPACE
