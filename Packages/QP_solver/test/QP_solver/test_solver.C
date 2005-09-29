// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
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
// file          : test/QP_solver/test_MPS.C
// package       : $CGAL_Package: _QP_solver $
//
// revision      : 0.2
// revision_date : 2000/08/21
//
// author(s)     : Kaspar Fischer (fischerk@inf.ethz.ch)
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test program for the QP solver
// ============================================================================

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <locale>

#include <CGAL/QP_solver/gmp_double.h>
#include <CGAL/QP_solver.h>
#include <CGAL/QP_full_exact_pricing.h>
#include <CGAL/QP_partial_exact_pricing.h>
#include <CGAL/QP_full_filtered_pricing.h>
#include <CGAL/QP_partial_filtered_pricing.h>
#include <CGAL/QP_solver/Double.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>

#include <CGAL/QP_solver/MPS.h> // should to into QP_solver.h (?)

// Note: The following #define's allow faster compilation for individual
// test cases.  For instance, to only compile code for a solver
// specialized for linear, nonsymmetric, has-equalities-only-and-full-rank,
// and standard form with internal integer arithmetic, pass
//
//   -DQP_L -DQP_NOT_S -DQP_R -DQP_F -DQP_INT
//
// to GNU g++.
#if (!defined(QP_L) && !defined(QP_S) && !defined(QP_R) && !defined(QP_F) && \
      !defined(QP_NOT_L) && !defined(QP_NOT_S) && !defined(QP_NOT_R) && \
      !defined(QP_NOTF))
#define QP_L     true
#define QP_S     true
#define QP_R     true
#define QP_F     true
#define QP_NOT_L true
#define QP_NOT_S true
#define QP_NOT_R true
#define QP_NOT_F true
#endif
#if (!defined(QP_INT) && !defined(QP_DOUBLE) && !defined(QP_RATIONAL))
#define QP_INT      true
#define QP_DOUBLE   true
#define QP_RATIONAL true
#endif

enum Pricing_strategy_type { FE, FF, PE, PF };
enum Input_type { Int_type, Double_type, Rational_type };
using CGAL::Tag_true;
using CGAL::Tag_false;
typedef std::map<std::string,int>::const_iterator Key_const_iterator;
typedef std::map<std::string,int>::iterator Key_iterator;
typedef std::pair<std::string,int> Arg;

// The following selector is used to derive the "cheap" number-type
// used during the pricing:
template<typename ET>
struct NT_selector {};

template<>
struct NT_selector<CGAL::Gmpq> {
  typedef CGAL::Gmpq NT;
};

template<>
struct NT_selector<CGAL::Double> {
  typedef double NT;
};

template<>
struct NT_selector<CGAL::Gmpz> {
  typedef CGAL::Gmpz NT;
};

template<typename T>
bool is_double(const T&)            { return false; }
bool is_double(const double&)       { return true; }
template<typename T>
bool is_int(const T&)               { return false; }
bool is_int(const int&)             { return true; }
template<typename T>
bool is_integer(const T&)           { return false; }
bool is_integer(const CGAL::Gmpz&)  { return true; }
template<typename T>
bool is_rational(const T&)          { return false; }
bool is_rational(const CGAL::Gmpq&) { return true; }

template <typename T>
T string_to(const std::string& s) {
  std::stringstream strm(s);
  T t;
  strm >> t;
  return t;
}

// Replaces the first occurrence of '%' in msg by replacement.
std::string replace1(const std::string& msg,const std::string& replacement)
{
  std::string result(msg);
  const std::string::size_type pos = result.find('%');
  CGAL_qpe_assertion(pos < result.size());
  result.replace(pos,1,replacement);
  return result;
}

void bailout(const char *msg)
{
  std::cout << "Error: " << msg << '.' << std::endl;
  exit(1);
}

void bailout1(const char *msg,const std::string& param)
{
  std::cout << "Error: " << replace1(msg,param) << '.' << std::endl;
  exit(1);
}

void usage()
{
  using std::cout;
  using std::endl;
  cout << "Normal usage:      ./test_solver < test_solver.cin\n"
       << "Alternative usage: ./test_solver verb strategy file [options]\n"
       << "\n"
       << "In the second form, verb is an integer (0-5), strategy is from\n"
       << "{fe,ff,pe,pf}, and file is a path to a MPS file. In addition,\n"
       << "you can specify any of the following additional options:\n"
       << "  +l       use the dedicated LP-solver on this instance\n"
       << "  +s       use the dedicated solver for instances whose\n"
       << "           D matrix is symmetric\n"
       << "  +r       use the dedicated solver for instances that only\n"
       << "           have equality constraints and whose coefficient\n"
       << "           matrix has full row rank\n"
       << "  +f       use the dedicated solver for instances in standard\n"
       << "           form (i.e., all variables have bounds [0,+inf]).\n"
       << "  int      assume the numbers in the MPS file are ints and use\n"
       << "           arbitrary-precision integers internally\n"
       << "  double   assume the numbers in the MPS file are doubles and\n"
       << "           use arbitrary-precision doubles internally\n"
       << "  rational assume the numbers in the MPS file are rationals\n"
       << "           use arbitrary-precision rationals internally\n"
       << "If any option is not specified, the program will test all\n"
       << "possible combinations (e.g., if '+s' is not specified, it will\n"
       << "check whether D is symmetric and run both, the dedicated solver\n"
       << "for symmetric and the non-dedicated solver on the instance).\n"
       << "You can also negate the options 'l', 's', 'r', or 'f' by\n"
       << "replacing the '-' by a '+'.\n";
  exit(1);
}

namespace Token { // A simple token reader with put-back facility.

  bool have_put_back_ = false;
  std::string put_back_token;

  std::string token(std::istream& in)
  {
    if (have_put_back_) {
      have_put_back_ = false;
      return put_back_token;
    }
    std::string t;
    in >> t;
    return t;
  }

  void put_token_back(const std::string& t)
  {
    have_put_back_ = true;
    put_back_token = t;
  }

} // namespace Token

// Clears options and reads from in the arguments as specified in usage().
// Returns true iff a new set of options could be parsed; if so, filename
// receives the name of the MPS file to solve.
bool parse_options(std::istream& in,std::map<std::string,int>& options,
		   std::string &filename)
{
  options.clear();
  std::ws(in); // read white-space
  if (in.eof())
    return false;

  // read verbosity:
  std::string t = Token::token(in);
  const int v = string_to<int>(t);
  if (v<0 || v>5)
    bailout("illegal verbosity");
  options.insert(Arg("Verbosity",v));

  // read strategy:
  std::string st = Token::token(in);
  Pricing_strategy_type type;
  if (st=="fe")
    type = FE;
  else if (st=="ff")
    type = FF;
  else if (st=="pe")
    type = PE;
  else if (st=="pf")
    type = PF;
  else
    bailout1("illegal pricing strategy '%'",st);
  options.insert(Arg("Strategy",static_cast<int>(type)));

  // read file name:
  in >> filename;
  if (filename.size() == 0)
    bailout("no filename specified");

  // read additional options:
  std::ws(in);
  bool eof = in.eof();
  t = Token::token(in);
  while (!eof && !isdigit(t[0])) {
    // process:
    bool good;
    if (t == "+l")
      good = options.insert(Arg("Is linear",1)).second;
    else if (t == "-l")
      good = options.insert(Arg("Is linear",0)).second;
    else if (t == "+s")
      good = options.insert(Arg("Is symmetric",1)).second;
    else if (t == "-s")
      good = options.insert(Arg("Is symmetric",0)).second;
    else if (t == "+r")
      good = options.insert(Arg("Has equalities and full rank",1)).second;
    else if (t == "-r")
      good = options.insert(Arg("Has equalities and full rank",0)).second;
    else if (t == "+f")
      good = options.insert(Arg("Is in standard form",1)).second;
    else if (t == "-f")
      good = options.insert(Arg("Is in standard form",0)).second;
    else {
      Input_type input;
      if (t == "int")
	input = Int_type;
      else if (t == "double")
	input = Double_type;
      else if (t == "rational")
	input = Rational_type;
      else
	bailout1("unknown input number type '%'",t);
      good = options.insert(Arg("Input type",static_cast<int>(input))).second;
    }

    // process duplicate errors:
    if (!good)
      bailout1("duplicate/contradictory option for '%'",t);

    // advance:
    ws(in);
    eof = in.eof();
    t = Token::token(in);
  }
  Token::put_token_back(t);

  // output:
  const int Width = 35;
  using std::cout;
  using std::left;
  using std::setw;
  using std::endl;
  cout << left << setw(Width) << "Processing:" << filename << endl;
  for (Key_const_iterator it = options.begin();
       it != options.end(); ++it) {
    cout << it->first << left << setw(Width - it->first.size()) <<  ":";
    if (it-> first == "Strategy")
      cout << st << endl;
    else
      cout << it->second << endl;
  }

  return true;
}

template<typename Traits>
CGAL::QP_pricing_strategy<Traits> *
  create_strategy(const std::map<std::string,int>& options)
{
  Key_const_iterator it = options.find("Strategy");
  CGAL::QP_pricing_strategy<Traits> *strat;
  typedef typename NT_selector<typename Traits::ET>::NT NT;
  switch (it->second) {
  case FE:
    strat = new CGAL::QP_full_exact_pricing<Traits>;
    break;
  case FF:
    strat = new CGAL::QP_full_filtered_pricing<Traits,NT>;
    break;
  case PE:
    strat = new CGAL::QP_partial_exact_pricing<Traits>;
    break;
  case PF:
    strat = new CGAL::QP_partial_filtered_pricing<Traits,NT>;
  }
  return strat;
}

template<typename Is_linear,
	 typename Is_symmetric,
	 typename Has_equalities_only_and_full_rank,
	 typename Is_in_standard_form,
	 typename IT,
	 typename ET>
bool process(const std::string& filename,
	     const std::map<std::string,int>& options)
{
  using std::cout;
  using std::endl;

  // extract verbosity:
  const int verbosity = options.find("Verbosity")->second;

  // read QP instance:
  std::ifstream in(filename.c_str());
  if (!in)
    bailout1("could not open file '%'",filename);
  typedef CGAL::QP_MPS_instance<IT,ET> QP_instance;
  QP_instance qp(in,true,verbosity);
  in.close();

  // check whether we should compute the rank of the coefficient
  // matrix:
  std::string comment = qp.comment();
  const bool dontComputeRank = 
    comment.find("dont-compute-row-rank") < comment.size();

  // check for the number-type:
  Input_type type;
  std::string number_type;
  if (comment.find("Number-type: floating-point") < comment.size()) {
    type = Double_type;
    number_type = "double";
  } else if (comment.find("Number-type: integer") < comment.size()) {
    type = Int_type;
    number_type = "Gmpz";
  } else if (comment.find("Number-type: rational") < comment.size()) {
    type = Rational_type;
    number_type = "Gmpq";
  } else
    bailout("file does not specify a number-type");
  if (type==Double_type && (is_int(IT()) || is_rational(IT())) ||
      type==Int_type && is_rational(IT()) ||
      type==Rational_type && (is_double(IT()) || is_int(IT())))
    return true;

  // construct a zero D matrix if needed:
  if (qp.is_linear() && !check_tag(Is_linear()))
    // Note: Revision 1.1 of this file uses qp's zero_D() routine and
    // a special traits class for the solver for this case.  But as
    // this more than doubles the compilation time, I removed it
    // again...
    qp.make_zero_D(); 

  // check which properties the loaded QP has, and break if they are
  // in contradiction to the routine's compile-time flags: 
  if (check_tag(Is_linear()) && !qp.is_linear() ||
      check_tag(Is_symmetric()) && !qp.is_symmetric() ||
      check_tag(Is_in_standard_form()) && !qp.is_in_standard_form() ||
      check_tag(Has_equalities_only_and_full_rank()) && dontComputeRank ||
      check_tag(Has_equalities_only_and_full_rank()) &&
      !qp.has_equalities_only_and_full_rank())
    return true;

  if (verbosity > 0)
    cout << "- Running a solver specialized for: "
	 << (check_tag(Is_linear())? "linear " : "")
	 << (check_tag(Is_symmetric())? "symmetric " : "")
	 << (check_tag(Is_in_standard_form())? "standard-form " : "")
	 << (check_tag(Has_equalities_only_and_full_rank())?
	     "has-equalities-only-and-full-rank " : "") 
	 << "file-IT=" << number_type << ' '
	 << "IT=" << (is_double(IT())? "double" :
		      (is_rational(IT())? "Gmpq" : "int")) << ' '
	 << "ET=" << (is_integer(ET())? "Gmpz" :
		      (is_rational(ET())? "Gmpq" : "Double"))
	 << endl;

  // check for format errors in MPS file:
  if (!qp.is_valid()) {
    cout << "Input is not a valid MPS file." << endl
	 << "Error: " << qp.error() << endl;
    return false;
  }
  if (verbosity > 3)
    cout << endl << qp;

  typedef CGAL::QP_solver_MPS_traits_d<QP_instance,
    Is_linear,Is_symmetric,Has_equalities_only_and_full_rank,
    Is_in_standard_form,IT,ET,
    typename QP_instance::D_iterator> Traits;

  // solve:
  CGAL::QP_pricing_strategy<Traits> *s = create_strategy<Traits>(options);
  CGAL::QP_solver<Traits> solver(qp.number_of_variables(),
				 qp.number_of_constraints(),
				 qp.A(),qp.b(),qp.c(),qp.D(),
				 qp.row_types(),
				 qp.fl(),qp.l(),qp.fu(),qp.u(),
				 s,verbosity);
  const bool is_valid = solver.is_valid();
  delete s;

  if (verbosity > 0 || !is_valid)
    cout << "  Solution is valid: " << is_valid << endl;
  return is_valid;
}

template<typename Is_linear,
	 typename Is_symmetric,
	 typename Has_equalities_only_and_full_rank,
	 typename Is_in_standard_form>
bool processType(const std::string& filename,
		 const std::map<std::string,int>& options)
{
  // look up value:
  Key_const_iterator it = options.find("Input type");
  const bool processOnlyOneValue = it != options.end();
  Input_type value;
  if (processOnlyOneValue)
    value = static_cast<Input_type>(it->second);

  // do only this particular value or all possibilities:
  bool success = true;
#ifdef QP_INT
  if (!processOnlyOneValue || value==Int_type)
    if (!process<Is_linear,Is_symmetric,
	Has_equalities_only_and_full_rank,Is_in_standard_form,
	int,CGAL::Gmpz>(filename,options))
      success = false;
#endif
#ifdef QP_DOUBLE
  if (!processOnlyOneValue || value==Double_type)
    if (!process<Is_linear,Is_symmetric,
	Has_equalities_only_and_full_rank,Is_in_standard_form,
	double,CGAL::Double>(filename,options))
      success = false;
#endif
#ifdef QP_RATIONAL
  if (!processOnlyOneValue || value==Rational_type)
    if (!process<Is_linear,Is_symmetric,
	Has_equalities_only_and_full_rank,Is_in_standard_form,
	CGAL::Gmpq,CGAL::Gmpq>(filename,options))
      success = false;
#endif
  return success;
}

template<typename Is_linear,
	 typename Is_symmetric,
	 typename Has_equalities_only_and_full_rank>
bool processFType(const std::string& filename,
		  const std::map<std::string,int>& options)
{
  Key_const_iterator it = options.find("Is in standard form");
  const bool processOnlyOneValue = it != options.end();
  bool value;
  if (processOnlyOneValue)
    value = it->second > 0;
  bool success = true;
#ifdef QP_F
  if (!processOnlyOneValue || value==true)
    if (!processType<Is_linear,Is_symmetric,
	Has_equalities_only_and_full_rank,Tag_true>(filename,options))
      success = false;
#endif
#if 0 // todo: once solver works for explicit bounds, we change this
      // back to: #ifdef QP_NOT_F
  if (!processOnlyOneValue || value==false)
    if (!processType<Is_linear,Is_symmetric,
	Has_equalities_only_and_full_rank,Tag_false>(filename,options))
      success = false;
#endif
  return success;  
}

template<typename Is_linear,
	 typename Is_symmetric>
bool processRFType(const std::string& filename,
		   const std::map<std::string,int>& options)
{
  Key_const_iterator it = options.find("Has equalities and full rank");
  const bool processOnlyOneValue = it != options.end();
  bool value;
  if (processOnlyOneValue)
    value = it->second > 0;
  bool success = true;
#ifdef QP_R
  if (!processOnlyOneValue || value==true)
    if (!processFType<Is_linear,Is_symmetric,Tag_true>(filename,options))
      success = false;
#endif
#ifdef QP_NOT_R
  if (!processOnlyOneValue || value==false)
    if (!processFType<Is_linear,Is_symmetric,Tag_false>(filename,options))
      success = false;
#endif
  return success;  
}

template<typename Is_linear>
bool processSRFType(const std::string& filename,
		    const std::map<std::string,int>& options)
{
  Key_const_iterator it = options.find("Is symmetric");
  const bool processOnlyOneValue = it != options.end();
  bool value;
  if (processOnlyOneValue)
    value = it->second > 0;
  bool success = true;
#ifdef QP_S
  if (!processOnlyOneValue || value==true)
    if (!processRFType<Is_linear,Tag_true>(filename,options))
      success = false;
#endif
#ifdef QP_NOT_S
  if (!processOnlyOneValue || value==false)
    if (!processRFType<Is_linear,Tag_false>(filename,options))
      success = false;
#endif
  return success;  
}

bool processLSRFType(const std::string& filename,
		     const std::map<std::string,int>& options)
{
  Key_const_iterator it = options.find("Is linear");
  const bool processOnlyOneValue = it != options.end();
  bool value;
  if (processOnlyOneValue)
    value = it->second > 0;
  bool success = true;
#ifdef QP_L
  if (!processOnlyOneValue || value==true)
    if (!processSRFType<Tag_true>(filename,options))
      success = false;
#endif
#ifdef QP_NOT_L
  if (!processOnlyOneValue || value==false)
    if (!processSRFType<Tag_false>(filename,options))
      success = false;
#endif
  return success;  
}

int main(const int ac,const char **av) {
  using std::cout;
  using std::endl;
  using CGAL::Tag_true;
  using CGAL::Tag_false;

  // determine in which mode we run:
  const bool readFromStdIn = ac == 1;
  std::stringstream args;
  if (!readFromStdIn) {
    if (ac < 4)
      usage();
    for (int i=1; i<ac; ++i) // convert av to a stream
      args << ' ' << av[i];
    cout << "Reading arguments from command line..." << endl;
  } else 
    cout << "Reading from standard in..." << endl;
  std::istream &in = (readFromStdIn? std::cin : args);

  // temporarily until upper bounding works:
  cout << "Warning: test_solver currently does not run" << endl
       << "the general solver for non-standard-form problems" << endl
       << "(because it is not yet fully functional)." << endl;

  // process input file(s):
  std::map<std::string,int> options;
  std::string filename;
  bool success = true;
  while (parse_options(in,options,filename))
    if (!processLSRFType(filename,options))
      success = false;

  // final output:
  if (!success) {
    cout << "Warning: some test cases were not handled correctly." << endl;
    return 1;
  } else {
    cout << "All test cases were successfully passed." << endl;
    return 0;
  }
}
