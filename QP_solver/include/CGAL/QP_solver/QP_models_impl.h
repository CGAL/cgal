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

#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>

CGAL_BEGIN_NAMESPACE

template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
QP_from_mps<IT_, Is_linear_,
		Sparse_D_,
		Sparse_A_>::
QP_from_mps(std::istream& in,bool use_CPLEX_convention,
		int verbosity)
  : has_linear_tag (check_tag(Is_linear_())),
    verbosity_(verbosity), from(in),
    is_format_okay_(false),
    is_linear_(true),
    use_CPLEX_convention(use_CPLEX_convention),
    it0(0),
    is_nonnegative_cached(false),
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

  // initialize matrix D (in the linear case, this is not needed,
  // since we return a Const_oneset_iterator that never accesses D
  if (!has_linear_tag)
    initialize_D(var_names.size(), Sparse_D());

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

template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
bool QP_from_mps<IT_, Is_linear_,
		Sparse_D_,
		Sparse_A_>::is_valid() const
{
  return is_format_okay_;
}

template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
const std::string& QP_from_mps<IT_, Is_linear_,
		Sparse_D_,
		Sparse_A_>::error()
{
  CGAL_qpe_assertion(!is_valid());
  return error_msg;
}

template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
const std::string& QP_from_mps<IT_, Is_linear_,
		Sparse_D_,
		Sparse_A_>::comment()
{
  return comment_;
}

template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
bool QP_from_mps<IT_, Is_linear_, Sparse_D_, Sparse_A_>
::is_symmetric(Tag_false /*sparse_D*/, unsigned int&i, unsigned int&j) const
{
  // only called if we have a qp, i.e. if D is initialized
  const unsigned int var_nr = var_names.size();
  for (i=0; i<var_nr; ++i)
    for (j=i+1; j<var_nr; ++j)
      if (D_[i][j] != D_[j][i])
	return false;
  return true;
}

template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
bool QP_from_mps<IT_, Is_linear_, Sparse_D_, Sparse_A_>
::is_symmetric(Tag_true /*sparse_D*/, unsigned int&i, unsigned int&j) const
{
  // only called if we have a qp, i.e. if D is initialized
  const unsigned int var_nr = var_names.size();
  for (i=0; i<var_nr; ++i) {
    typename Map::const_iterator b = D_[i].begin();
    typename Map::const_iterator e = D_[i].end();
    for (; b != e; ++b) {
      j = b->first;
      // check entry (j, i)
      IT expected_val = b->second;
      IT found_val = it0;
      typename Map::const_iterator it = D_[j].find(i);
      if (it !=  D_[j].end()) found_val = it->second;
      if (expected_val != found_val) return false;
    }
  }
  return true;
}
  

template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
bool QP_from_mps<IT_, Is_linear_,
		Sparse_D_,
		Sparse_A_>::name_section()
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

template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
bool QP_from_mps<IT_, Is_linear_,
		Sparse_D_,
		Sparse_A_>::rows_section()
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

template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
bool QP_from_mps<IT_, Is_linear_,
		Sparse_D_,
		Sparse_A_>::columns_section()
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
      add_column (Sparse_A_());
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
	set_entry_in_A (var_index, row_name->second, val, Sparse_A_());
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

template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
bool QP_from_mps<IT_, Is_linear_,
		Sparse_D_,
		Sparse_A_>::rhs_section()
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

template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
bool QP_from_mps<IT_, Is_linear_,
		Sparse_D_,
		Sparse_A_>::ranges_section()
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
	  int new_index = b_.size();
	  for (unsigned int j=0; j<var_names.size(); ++j) {
	    IT val = get_entry_in_A (j, index, Sparse_A_());
	    add_entry_in_A (j, new_index, val, Sparse_A_());  
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
	 


template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
bool QP_from_mps<IT_, Is_linear_,
		Sparse_D_,
		Sparse_A_>::bounds_section()
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

template<typename IT_, typename Is_linear_,
	 typename Sparse_D_,
	 typename Sparse_A_>
bool QP_from_mps<IT_, Is_linear_,
		Sparse_D_,
		Sparse_A_>::qmatrix_section()
{
  std::string t = token();
  if (t!="QMATRIX" && t!="DMATRIX" && t!="QUADOBJ") { // optional
    put_token_back(t);
    return true;
  }   

  // remember section name:
  D_section = t;
  const bool multiply_by_two = t=="DMATRIX";
  const bool only_get_lower_part = t =="QUADOBJ";

  t = token();
  std::string bound_id;
  while (t != "ENDATA") {
    // find first variable name;
    const Index_map::const_iterator var1_name = var_names.find(t);
    if (var1_name == var_names.end()) // unknown variable?
      return err2("unknown first variable '%' in '%' section", t, D_section);
    const unsigned int var1_index = var1_name->second;;
    //std::cout << "qvar1 " << t << std::endl;
      
    // find second variable name;
    t = token();
    const Index_map::const_iterator var2_name = var_names.find(t);
    if (var2_name == var_names.end()) // unknown variable?
      return err2("unknown second variable '%' in '%' section",t, D_section);
    const unsigned int var2_index = var2_name->second;;
    //std::cout << "qvar2 " << t << std::endl;
      
    // read value:
    IT val;
    if (!number(val))
      return err1("expected number after '%' in section QMATRIX",t);

    // multiply by two if approriate:
    if (multiply_by_two)
      val *= IT(2);

    // mark problem as nonlinear if value is nonzero, and
    // bail out if we are supposed to read a linear program
    if (is_linear_ && !CGAL::is_zero(val)) {
      is_linear_ = false;
      if (has_linear_tag)
	return err1
	  ("nonzero value in '%' section while reading linear program", 
	   D_section);
    }

    // set entry in D:
    set_entry_in_D(var1_index, var2_index, val, Sparse_D());
    if (only_get_lower_part)
      // duplicate entry if not on diagonal
      if (var1_index != var2_index) 
	set_entry_in_D(var2_index, var1_index, val, Sparse_D());

    // read next token:
    t = token();
  }
  put_token_back(t);

  // now check symmetry
  unsigned int bad_i, bad_j;
  if (!is_symmetric(Sparse_D(), bad_i, bad_j)) {
    std::string bad_i_name = var_by_index[bad_i];
    std::string bad_j_name = var_by_index[bad_j];
    return 
      err3("nonsymmetric '%' section for variables '%' and '%'",
	   D_section, bad_i_name, bad_j_name);
  }
	   
  return true;
}

CGAL_END_NAMESPACE
