// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <locale>
#include <map>
#include <cstdlib>
#include <CGAL/config.h>

#ifndef CGAL_USE_GMP
#include <CGAL/MP_Float.h>   // with normalization switched on
#else
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpzf.h>
#endif
#include <CGAL/Exact_rational.h>

#include <CGAL/QP_solver/QP_solver.h>
#include <CGAL/QP_solver/QP_full_exact_pricing.h>
#include <CGAL/QP_solver/QP_exact_bland_pricing.h>
#include <CGAL/QP_solver/QP_partial_exact_pricing.h>
#include <CGAL/QP_solver/QP_full_filtered_pricing.h>
#include <CGAL/QP_solver/QP_partial_filtered_pricing.h>
#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>

// Note: The following #define's allow faster compilation for individual
// test cases.  For instance, to only compile code for a solver
// specialized for linear, nonsymmetric, has-equalities-only-and-full-rank,
// and standard form with internal integer arithmetic, pass
//
//   -DQP_L -DQP_NOT_S -DQP_R -DQP_F -DQP_INT
//
// to GNU g++.  (Or even simpler, work with test_MPS.cpp.)
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

enum Input_type { Int_type, Double_type, Rational_type };
using CGAL::Tag_true;
using CGAL::Tag_false;
typedef std::map<std::string,int>::const_iterator Key_const_iterator;
typedef std::map<std::string,int>::iterator Key_iterator;
typedef std::pair<std::string,int> Arg;

#ifdef CGAL_USE_GMP
typedef CGAL::Gmpzf Float;
typedef CGAL::Gmpz Integer;
#else
typedef CGAL::MP_Float Float;
typedef CGAL::MP_Float Integer;
//typedef CGAL::Quotient<CGAL::MP_Float> Rational;  // too slow
#endif
typedef CGAL::Exact_rational Rational;

template<typename T>
bool is_double(const T&)            { return false; }
bool is_double(const double&)       { return true; }
template<typename T>
bool is_int(const T&)               { return false; }
bool is_int(const int&)             { return true; }
template<typename T>
bool is_Integer(const T&)           { return false; }
bool is_Integer(const Integer&)     { return true; }
template<typename T>
bool is_rational(const T&)          { return false; }
bool is_rational(const Rational&)   { return true; }

int string_to_int(const std::string& s) {
  std::stringstream strm(s);
  int t;
  strm >> t;
  return t;
}

// Replaces the first occurrence of '%' in msg by replacement.
std::string replace1(const std::string& msg,const std::string& replacement)
{
  std::string result(msg);
  const std::string::size_type pos = result.find('%');
  assert(pos < result.size());
  result.replace(pos,1,replacement);
  return result;
}

void bailout(const char *msg)
{
  std::cout << "Error: " << msg << '.' << std::endl;
  std::exit(1);
}

void bailout1(const char *msg,const std::string& param)
{
  std::cout << "Error: " << replace1(msg,param) << '.' << std::endl;
  std::exit(1);
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
       << "  +f       use the dedicated solver for instances in standard\n"
       << "           form (i.e., all variables have bounds [0,+inf]).\n"
       << "  int      assume the numbers in the MPS file are ints and use\n"
       << "           arbitrary-precision integers internally\n"
       << "  double   assume the numbers in the MPS file are doubles and\n"
       << "           use arbitrary-precision doubles internally\n"
       << "  rational assume the numbers in the MPS file are rationals\n"
       << "           use arbitrary-precision rationals internally\n"
       << "If any option is not specified, the program will test all\n"
       << "possible combinations (e.g., if '+l' is not specified, it will\n"
       << "check whether the problem is linear and run both, the dedicated\n"
       << "solver for the linear case and the non-dedicated solver on the\n"
       << " instance).\n"
       << "You can also negate the options 'l'or 'f' by\n"
       << "replacing the '-' by a '+'.\n";
  std::exit(1);
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
  const int v = string_to_int(t);
  if (v<0 || v>5)
    bailout("illegal verbosity");
  options.insert(Arg("Verbosity",v));

  // read strategy:
  std::string st = Token::token(in);
  CGAL::Quadratic_program_pricing_strategy type = CGAL::QP_CHOOSE_DEFAULT;
  if (st=="fe")
    type = CGAL::QP_DANTZIG;
  else if (st=="eb")
    type = CGAL::QP_BLAND;
  else if (st=="ff")
    type = CGAL::QP_FILTERED_DANTZIG;
  else if (st=="pe")
    type = CGAL::QP_PARTIAL_DANTZIG;
  else if (st=="pf")
    type = CGAL::QP_PARTIAL_FILTERED_DANTZIG;
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
    else if (t == "+f")
      good = options.insert(Arg("Is in standard form",1)).second;
    else if (t == "-f")
      good = options.insert(Arg("Is in standard form",0)).second;
    else {
      Input_type input = Int_type; // to kill warnings
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
  const int Width = 15;
  using std::cout;
  using std::left;
  using std::setw;
  using std::endl;
  cout << endl << left << setw(Width) << "Processing:" << filename << endl;
  for (Key_const_iterator it = options.begin();
       it != options.end(); ++it) {
    cout << " " << it->first << left << setw(Width-it->first.size()-1) <<  ":";
    if (it-> first == "Strategy")
      cout << st << endl;
    else
      cout << it->second << endl;
  }

  return true;
}

template<typename Is_linear,
         typename Is_nonnegative,
         typename IT,
         typename ET>
bool process(const std::string& filename,
             const std::map<std::string,int>& options)
{
  using std::cout;
  using std::endl;

  // extract verbosity:
  const int verbosity = options.find("Verbosity")->second;

  // read QP instance
  std::ifstream in(filename.c_str());
  if (!in)
    bailout1("could not open file '%'",filename);
  typedef CGAL::Quadratic_program_from_mps <IT> QP_instance;
  QP_instance qp(in);
  in.close();

  std::string comment = qp.get_comment();

  // check for the number-type:
  Input_type type;
  std::string number_type;
  if (comment.find("Number-type: integer") < comment.size()) {
    type = Int_type;
    number_type = "Gmpz";
  } else if (comment.find("Number-type: rational") < comment.size()) {
    type = Rational_type;
    number_type = "Gmpq";
  } else { // "Number-type: floating-point" or none specifier:
    type = Double_type;
    number_type = "double";
  }
  // now, some combinations of IT and the file's number-type are
  // incomaptible:
  //    file's input type  | input type IT to be used in parsing the file
  //    -----------------------------------------------------------------
  //    double             | int  (can't convert double to it)
  //    int                |
  //    rational           | int + double (can't convert rational to them)
  if ((type==Double_type && is_int(IT()) ) ||
      (type==Rational_type && (is_double(IT()) || is_int(IT()))))
    return true;

  // also, Is_linear/Is_nonnegative tags are incompatible with
  // nonlinear/nonnegative programs
  if ((check_tag(Is_linear()) && !qp.is_linear()) ||
      (check_tag(Is_nonnegative()) && !qp.is_nonnegative()))
    return true;

  // finally, filtered pricing strategies are only to be used if the input
  // type is double
  Key_const_iterator it = options.find("Strategy");
  if (!is_double(IT()) &&
      (it->second == CGAL::QP_FILTERED_DANTZIG ||
       it->second == CGAL::QP_PARTIAL_FILTERED_DANTZIG))
    return true;

  if (verbosity > 0)
    cout << "- Running a solver specialized for: "
         << (check_tag(Is_linear())? "linear " : "")
         << (check_tag(Is_nonnegative())? "standard-form " : "")
         << "file-IT=" << number_type << ' '
         << "IT=" << (is_double(IT())? "double" :
                      (is_rational(IT())? "Gmpq" : "int")) << ' '
         << "ET=" << (is_Integer(ET())? "Gmpz" :
                      (is_rational(ET())? "Gmpq" : "Double"))
         << endl;

  // check for format errors in MPS file:
  if (!qp.is_valid()) {
    cout << "Input is not a valid MPS file." << endl
         << "Error: " << qp.get_error() << endl;
    return false;
  }

  // check name <-> index mapping
  for (int j=0; j<qp.get_n(); ++j)
    if (j != qp.variable_index_by_name (qp.variable_name_by_index (j))) {
      cout << "incorrect name <-> index mapping (variables)" << endl;
      return false;
    }
  for (int i=0; i<qp.get_m(); ++i)
    if (i != qp.constraint_index_by_name (qp.constraint_name_by_index (i))) {
      cout << "incorrect name <-> index mapping (constraints)" << endl;
      return false;
    }

  // print program (using QMATRIX format), read it back in and check
  // whether it still agrees with the original program
  std::stringstream inout;
  // if we have doubles, adjust precision to accommodate high-precision doubles
  if (is_double(IT())) inout << std::setprecision (12);
  CGAL::QP_functions_detail::print_program
    (inout, qp, std::string("test_io_mps"),
     Is_linear(),Is_nonnegative());

  // test readers
  CGAL::Quadratic_program_from_mps<IT> qp2(inout);
  assert (qp2.is_valid());
  if (!CGAL::QP_functions_detail::are_equal_qp (qp, qp2)) {
    cout << "Warning: MPS reader (QP) and MPS writer disagree.\n" << endl;
  }

  // now copy from the iterators, check for equality
  CGAL::Quadratic_program<IT>
    qp4 (qp.get_n(), qp.get_m(), qp.get_a(), qp.get_b(), qp.get_r(),
         qp.get_fl(), qp.get_l(), qp.get_fu(), qp.get_u(),
         qp.get_d(), qp.get_c(), qp.get_c0());
  if (!CGAL::QP_functions_detail::are_equal_qp (qp, qp4)) {
    cout << "Program not correctly copied.\n" << endl;
    return false;
  }
  // test consistency of types
  if (qp.is_linear() != qp4.is_linear()) {
    cout << "Program types inconsistent (linearity): "
         << qp.is_linear() << " vs. " << qp4.is_linear()
         << "\n"<< endl;
    return false;
  }
  if (qp.is_nonnegative() != qp4.is_nonnegative()) {
    cout << "Program types inconsistent (nonnegativity): "
         << qp.is_nonnegative() << " vs. " << qp4.is_nonnegative()
         << "\n"<< endl;
    return false;
  }


  // now comes the actual, output-generating run
  // make options
  CGAL::Quadratic_program_options solver_options;
  solver_options.set_verbosity(verbosity);
  solver_options.set_pricing_strategy
    (static_cast<CGAL::Quadratic_program_pricing_strategy>
     (options.find("Strategy")->second));
  solver_options.set_auto_validation(true);

  CGAL::Quadratic_program_solution<ET> solution =
    CGAL::QP_functions_detail::solve_program
    (qp, ET(0), Is_linear(), Is_nonnegative(), solver_options);
  // output solution + number of iterations
  cout << CGAL::to_double(solution.objective_value()) << "("
       << solution.number_of_iterations() << ") ";
  bool is_valid = solution.is_valid();

  // the last step: solve poblem from copied QP, if full exact pricing
  // and general form
  if (options.find("Strategy")->second == CGAL::QP_DANTZIG &&
      !check_tag(Is_linear()) &&
      !check_tag(Is_nonnegative()))
    {
    // general form
    typedef CGAL::Quadratic_program<IT>
      LocalQP;
    CGAL::Quadratic_program_options local_options;
    local_options.set_verbosity(0);
    local_options.set_pricing_strategy(CGAL::QP_DANTZIG);
    local_options.set_auto_validation(true);
    LocalQP qplocal (qp.get_n(), qp.get_m(), qp.get_a(), qp.get_b(),
                     qp.get_r(),
                     qp.get_fl(), qp.get_l(), qp.get_fu(), qp.get_u(),
                     qp.get_d(), qp.get_c(), qp.get_c0());
    CGAL::solve_quadratic_program (qplocal, ET(), local_options);
    std::cout << "(c) ";
  }

  if (verbosity > 0 || !is_valid)
    cout << "  Solution is valid: " << is_valid << endl;
  return is_valid;
}

template<typename Is_linear,
         typename Is_nonnegative>
bool processType(const std::string& filename,
                 const std::map<std::string,int>& options)
{
  // look up value:
  Key_const_iterator it = options.find("Input type");
  const bool processOnlyOneValue = it != options.end();
  Input_type value = Int_type;
  if (processOnlyOneValue)
    value = static_cast<Input_type>(it->second);

  // do only this particular value or all possibilities:
  bool success = true;

#ifdef QP_INT
  if (!processOnlyOneValue || value==Int_type)
    if (!process<Is_linear,Is_nonnegative,
        int,Integer>(filename,options))
      success = false;
#endif
#ifdef QP_DOUBLE
  if (!processOnlyOneValue || value==Double_type)
    if (!process<Is_linear,Is_nonnegative,
        double,Float>(filename,options))
      success = false;
#endif
#ifdef QP_RATIONAL
  if (!processOnlyOneValue || value==Rational_type)
    if (!process<Is_linear,Is_nonnegative,
        Rational,Rational>(filename,options))
      success = false;
#endif

  return success;
}

template<typename Is_linear>
bool processFType(const std::string& filename,
                  const std::map<std::string,int>& options)
{
  Key_const_iterator it = options.find("Is in standard form");
  const bool processOnlyOneValue = it != options.end();
  bool value = false;
  if (processOnlyOneValue)
    value = it->second > 0;
  bool success = true;
#ifdef QP_F
  if (!processOnlyOneValue || value==true)
    if (!processType<Is_linear,Tag_true>(filename,options))
      success = false;
#endif
#ifdef QP_NOT_F
  if (!processOnlyOneValue || value==false)
    if (!processType<Is_linear,Tag_false>(filename,options))
      success = false;
#endif
  return success;
}

bool processLFType(const std::string& filename,
                   const std::map<std::string,int>& options)
{
  std::cout << " Solution:   ";
  Key_const_iterator it = options.find("Is linear");
  const bool processOnlyOneValue = it != options.end();
  bool value = false;
  if (processOnlyOneValue)
    value = it->second > 0;
  bool success = true;
#ifdef QP_L
  if (!processOnlyOneValue || value==true)
    if (!processFType<Tag_true>(filename,options))
      success = false;
#endif
#ifdef QP_NOT_L
  if (!processOnlyOneValue || value==false)
    if (!processFType<Tag_false>(filename,options))
      success = false;
#endif
  return success;
}

int main(const int ac,const char **av) {
  using std::cout;
  using std::endl;
  using CGAL::Tag_true;
  using CGAL::Tag_false;

  std::cerr.precision(17);
  std::cout.precision(17);

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
  std::istream in(readFromStdIn? std::cin.rdbuf() : args.rdbuf());

  // process input file(s):
  std::map<std::string,int> options;
  std::string filename;
  bool success = true;
  while (parse_options(in,options,filename))
    if (!processLFType(filename,options))
      success = false;

  // final output:
  cout << endl;
  if (!success) {
    cout << "Warning: some test cases were not handled correctly." << endl;
    return 1;
  } else {
    cout << "All test cases were successfully passed." << endl;
    return 0;
  }
}
