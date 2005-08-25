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

#include <CGAL/QP_partial_filtered_pricing.h>
#include <CGAL/QP_partial_exact_pricing.h>

CGAL_BEGIN_NAMESPACE

template<class Traits>
QP_MPS_instance<Traits>::QP_MPS_instance(std::istream& in,
					 bool use_CPLEX_convention,
					 QP_pricing_strategy<Traits> *strategy)
  : from(in),
    is_format_okay_(false),
    is_linear_(false),
    use_CPLEX_convention(use_CPLEX_convention),
    var_nr(0),
    constr_nr(0),
    strategy(strategy),
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

  std::cout << "MPS stream successfully read." << std::endl;
  is_format_okay_ = true;
}

template<class Traits>
QP_MPS_instance<Traits>::~QP_MPS_instance()
{
  if (strategy != 0)
    delete strategy;
}

template<class Traits>
bool QP_MPS_instance<Traits>::is_valid()
{
  if (!is_format_okay_)
    return false;

  // additional safety checks:
  // (Note: we do not check Has_equalities_only_and_full_rank currently,
  // as it is too expensive -- maybe add it for debugging purposes?)
  if (check_tag(Is_symmetric()) && !is_symmetric())
    return err("loaded instance does not have a symmetric D matrix but Is_symmetric is Tag_true");
  if (check_tag(Is_linear()) && !is_linear())
    return err("loaded instance is not an LP (D is nonzero) but Is_linear is Tag_true");
  if (check_tag(Is_in_standard_form()) && !is_in_standard_form())
    return err("loaded instance is not in standard form but Is_in_standard_form is Tag_true");
  return true;
}

template<class Traits>
const std::string QP_MPS_instance<Traits>::error()
{
  CGAL_qpe_assertion(!is_valid());
  return error_msg;
}

template<class Traits>
bool QP_MPS_instance<Traits>::is_symmetric()
{
  for (unsigned int i=0; i<var_nr; ++i)
    for (unsigned int j=i+1; j<var_nr; ++j)
      if (D_[i][j] != D_[j][i])
	return false;
  return true;
}
  
template<class Traits>
unsigned int QP_MPS_instance<Traits>::number_of_variables()
{
  CGAL_qpe_assertion(is_valid());
  return var_names.size();
}

template<class Traits>
unsigned int  QP_MPS_instance<Traits>::number_of_constraints()
{
  CGAL_qpe_assertion(is_valid());
  return row_names.size();
}

template<class Traits>
typename QP_MPS_instance<Traits>::A_iterator 
QP_MPS_instance<Traits>::A()
{
  CGAL_qpe_assertion(is_valid());
  return Vector_iterator(A_.begin(),Beginner());
}

template<class Traits>
typename QP_MPS_instance<Traits>::B_iterator
QP_MPS_instance<Traits>::b()
{
  CGAL_qpe_assertion(is_valid());
  return b_.begin();
}

template<class Traits>
typename QP_MPS_instance<Traits>::C_iterator
QP_MPS_instance<Traits>::c()
{
  CGAL_qpe_assertion(is_valid());
  return c_.begin();
}

template<class Traits>
typename QP_MPS_instance<Traits>::D_iterator
QP_MPS_instance<Traits>::D()
{
  CGAL_qpe_assertion(is_valid());
  return Vector_iterator(D_.begin(),Beginner());
}

template<class Traits>
typename QP_MPS_instance<Traits>::Row_type_iterator
QP_MPS_instance<Traits>::row_types()
{
  CGAL_qpe_assertion(is_valid());
  return row_types_.begin();
}

template<class Traits>
bool QP_MPS_instance<Traits>::name_section()
{
  const std::string t = token();
  if (t != "NAME")
    return err("expected 'NAME'");
  name = token();
  return true;
}

template<class Traits>
bool QP_MPS_instance<Traits>::rows_section()
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
      // register name of >=, <=, or = constraint:
      const unsigned int index = row_types_.size();
      row_types_.push_back(type == 'G'? Traits::GREATER_EQUAL :
			   (type == 'E'? Traits::EQUAL : Traits::LESS_EQUAL));
      if (row_names.find(t) != row_names.end())
	return err1("duplicate row name '%' in section ROWS",t);
      row_names.insert(String_int_pair(t,index));
      b_.push_back(IT(0));
      break;
    default:
      return err1(
		  "expected 'N', 'L', 'E', or 'G' in ROWS section but found '%'",
		  symbol);
    }
    t = token();
  }
  put_token_back(t);

  // remember number of constraints:
  constr_nr = row_names.size();
  return true;
}

template<class Traits>
bool QP_MPS_instance<Traits>::columns_section()
{
  std::string t = token();
  if (t != "COLUMNS")
    return err1("expected 'COLUMNS' but found '%'",t);

  t = token();
  while (t != "RHS") {
    // find variable name:
    unsigned int var_index;
    const Index_map::const_iterator var_name = var_names.find(t);
    if (var_name == var_names.end()) { // new variable?
      var_index = var_names.size();
      var_names.insert(String_int_pair(t,var_index));
      A_.push_back(Vector(constr_nr,IT(0)));
      c_.push_back(IT(0));
      fl_.push_back(true);  // default lower bound is finite...
      l_.push_back(IT(0));  // ...namely zero
      fu_.push_back(false); // default upper bound is infinite
      u_.push_back(IT());   // (dummy value) 
    } else // variable that is already known?
      var_index = var_name->second;
    //std::cout << "var is " << t << std::endl;
      
    bool doItAgain = true;
    for (int i=0; doItAgain; ++i) {
      // read row identifier:
      t = token();
      //std::cout << "row is " << t << std::endl;

      // read number:
      IT val;
      if (!number(val))
	return err1("number expected after row identifier '%' in this COLUMNS record",t);
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

template<class Traits>
bool QP_MPS_instance<Traits>::rhs_section()
{
  std::string t = token();
  if (t != "RHS")
    return err1("expected 'RHS' but found '%'",t);

  t = token();
  std::string rhs_id;
  while (t != "RANGES" && t != "BOUNDS" && t != "ENDATA") {
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

template<class Traits>
bool QP_MPS_instance<Traits>::bounds_section()
{
  std::string t = token();
  if (t != "BOUNDS") { // (Note: BOUNDS section is optional.)
    put_token_back(t);
    return true;
  }

  t = token();
  std::string bound_id;
  while (t != "QMATRIX" && t != "ENDATA") {
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
      type = FX;
    else if (t=="MI")
      type = MI;
    else if (t=="PL")
      type = PL;
    else    
      return
	err1("expected 'LO', 'UP', 'FX', 'FR', 'MI', or 'PL' here but found '%'",t);

    // find bound identifier:
    t = token();
    if (bound_id.size() == 0) // first time we see a bound identifier?
      bound_id = t;
    else if (t != bound_id) {
      warn1("ignoring all bounds for bound vector '%'",t);
      warn1("(only bounds for bound vector '%' are accepted)",bound_id);
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
	return err1("expected number after '%' in LO/UP/FX bound",t);

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

template<class Traits>
bool QP_MPS_instance<Traits>::qmatrix_section()
{
  std::string t = token();
  if (t != "QMATRIX") { // (Note: QMATRIX section is optional.)
    put_token_back(t);
    is_linear_ = true;
    return true;
  }

  // initialize matrix D:
  initialize_D(var_names.size(),Use_sparse_representation_for_D());

  t = token();
  std::string bound_id;
  while (t != "ENDDATA") {
    // find first variable name;
    const Index_map::const_iterator var1_name = var_names.find(t);
    if (var1_name == var_names.end()) // unknown variable?
      return err1("unknown first variable '%' in QMATRIX section",t);
    const unsigned int var1_index = var1_name->second;;
      
    // find second variable name;
    t = token();
    const Index_map::const_iterator var2_name = var_names.find(t);
    if (var2_name == var_names.end()) // unknown variable?
      return err1("unknown second variable '%' in QMATRIX section",t);
    const unsigned int var2_index = var2_name->second;;
      
    // read value of bound, if appropriate:
    IT val;
    if (!number(val))
      return err1("expected number after '%' in section QMATRIX",t);

    // set entry in D:
    set_entry_in_D(var1_index,var2_index,val,
		   Use_sparse_representation_for_D());

    // read next token:
    t = token();
  }
  put_token_back(t);

  return true;
}

template<class Traits>
std::ostream& operator<<(std::ostream& o,QP_MPS_instance<Traits>& qp)
{
  return o;
}

CGAL_END_NAMESPACE
