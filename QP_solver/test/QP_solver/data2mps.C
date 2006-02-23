// Synopsis: converts a .data-file to a MPS-file
//
// Usage: ./data2mps < data-file
//
// Author: Kaspar Fischer <fischerk@inf.ethz.ch>
//
// $Revision: $
// $Id: $

#include <iostream>
#include <sstream>
#include <string>

#include <boost/any.hpp>

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

class Data_reader {
public: // types:
  typedef boost::any ANT;
  typedef CGAL::Gmpq Rational;
  enum Row_type { LESS_EQUAL = -1, EQUAL, GREATER_EQUAL};

private: // types:
  enum Input_number_type { INTEGER, RATIONAL, DOUBLE, YET_UNKOWN_NUMBER_TYPE };

  typedef std::vector<ANT>                   Vector;
  typedef std::vector<Vector>                Matrix;
  typedef CGAL::QP_MPS_detail::Begin<Vector> Beginner;
  typedef CGAL::Join_input_iterator_1<Matrix::const_iterator,
			      Beginner >     Vector_iterator;
  typedef Vector::const_iterator             Entry_iterator;
  typedef std::vector<bool>                  F_vector;
  typedef F_vector::const_iterator           F_vector_iterator;
  typedef std::vector<Row_type>              Row_type_vector;

public: // types:

  // iterators over the input matrices and vectors:
  typedef Vector_iterator   A_iterator;
  typedef Entry_iterator    B_iterator;
  typedef Entry_iterator    C_iterator;
  typedef Vector_iterator   D_iterator;
  typedef CGAL::Const_oneset_iterator< CGAL::Const_oneset_iterator<ANT> >
                            Zero_D_iterator;
  typedef F_vector_iterator FU_iterator;
  typedef F_vector_iterator FL_iterator;
  typedef Entry_iterator    U_iterator;
  typedef Entry_iterator    L_iterator;

  typedef Row_type_vector::const_iterator Row_type_iterator;

private: // members:
  std::istream& in;               // input stream containing .data format
  bool good_;                     // whether or not we have successfully
				  // parsed the input
  std::string error_message_;     // error message, see error_message()

  // arguments read from the HEADER section:
  Input_number_type nt;           // what kind of numbers we are going to read
  ANT any_zero;                   // the number 0
  int n, m;                       // number of variables and constraints
  bool default_fl, default_fu;    // "is lower/upper bound finite" flags for
                                  // default bounds
  ANT default_l, default_u;       // default lower/upper bound

  // data:
  Matrix A_, D_;
  Vector b_, c_;
  F_vector fl_, fu_;
  Vector l_, u_;
  Row_type_vector row_type_;

public:
  Data_reader(std::istream& in) : in(in), good_(false),
				  nt(YET_UNKOWN_NUMBER_TYPE),
				  n(-1), m(-1)
  {
    if (!read_initial_comments()) {
      err("Unexpected EOF before HEADER");
      return;
    }

    if (!read_HEADER())
      return;

    // allocate storage:
    A_.resize(n, Vector(m, any_zero));
    D_.resize(n, Vector(n, any_zero));
    fl_.resize(n, default_fl);
    fu_.resize(n, default_fu);
    l_.resize(n, default_l);
    u_.resize(n, default_u);
    b_.resize(n);
    c_.resize(n);
    row_type_.resize(n);

    if (!read_CONSTRAINTS())
      return;

    std::string tk;
    if (!read_D_MATRIX(tk))
      return;

    if (!read_C_VECTOR_AND_BOUNDS(tk))
      return;

    whitespace();
    if (token() != "END") {
      err("Expected END after C-VECTOR-AND-BOUNDS section");
      return;
    }

    // having survived everything till here, we mark the input as valid:
    good_ = true;
  }

  bool good() const
  {
    return good_;
  }

  std::string& error_message()
  {
    CGAL_qpe_assertion(!good());
    return error_message_;    
  }

private: // error handling helpers:
  // Note: alternatively, I could have used exceptions to handle errors...

   template <typename T>
   std::string tostr(const T& t) {
     std::stringstream strm;
     strm << t;
     return strm.str();
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

  bool err(const char* msg)
  {
    error_message_ = msg;
    return false;
  }

  bool err1(const char* msg,
	    const std::string& parameter1)
  {
    error_message_ = replace1(msg,parameter1);
    return false;
  }

  bool err2(const char* msg,
	    const std::string& parameter1,
	    const std::string& parameter2)
  {
    error_message_ = replace1(replace1(msg,parameter1),parameter2);
    return false;
  }

  bool err3(const char* msg,
	    const std::string& parameter1,
	    const std::string& parameter2,
	    const std::string& parameter3)
  {
    error_message_ = replace1(replace1(replace1(msg,parameter1),parameter2),
			      parameter3);
    return false;
  }

  void warn(const std::string& msg)
  {
    std::cerr << "Parser warning: " << msg << '.' << std::endl;
  }

  void warn1(const std::string& msg,const std::string& parameter1)
  {
    warn(replace1(msg,parameter1));
  }

private: // parsing helpers:

  // Eats whitespace and returns true iff the end of the file has not yet been
  // reached.
  bool whitespace()
  {
    std::ws(in);
    return !in.eof();
  }

  // Returns the next token, i.e., string of consecutive non-whitespace
  // characters.
  std::string token() {
    whitespace();
    std::string token;
    char c;
    while (in.get(c)) {
      if (isspace(c)) {
	in.putback(c);
	break;
      }
      token.push_back(c);
    }

    // std::cout << "[token: '" << token << "']\n";
    return token;
  }

  bool read_number(ANT& val)
  {
    if (nt == INTEGER) {
      int v;
      in >> v;
      val = v;
    } else if (nt == RATIONAL) {
      Rational v;
      in >> v;
      val = v;
    } else  if (nt == DOUBLE) {
      double v;
      in >> v;
      val = v;
    } else
      CGAL_qpe_assertion(false);
    return in.good();
  }

private: // parsing routines:

  // Reads the comments preceeding HEADER, returning false in case of an error.
  // Note: this method does not touch error_message_.
  bool read_initial_comments()
  {
    do {
      // obtain first (real) charactor of the line:
      whitespace();
      int first_char = in.peek();
      if (first_char != '#')
	return in.good();

      // skip this line:
      while (in.good() && in.get() != '\n')
	;
    } while (in.good());
    return false;
  }

  // Reads the header, returning false in case of an error.
  // Note: this method sets error_message_ appropriately in case of errors.
  bool read_HEADER() {
    // safety check:
    std::string tk = token();
    if (tk != "HEADER")
      return err1("Expected HEADER after initial comments but found '%'", tk);

    // read header tags: 
    tk = token();
    while (tk != "CONSTRAINTS") {
      if (tk == "number-type:") {
	std::string nts = token();
	if (nts == "integer") {
	  nt = INTEGER;
	  any_zero = int(0);
	} else if (nts == "rational") {
	  nt = RATIONAL;
	  any_zero = Rational(0);
	} else if (nts == "double") { 
	  nt = DOUBLE;
	  any_zero = double(0);
	} else
	  return err1("Unknown type '%' in 'number-type:' clause in "
		      "HEADER section", nts);
	
	// set default bound:
	default_fl = true;
	default_fu = false;
	default_l  = any_zero;

      } else if (tk == "number-of-variables:") {
	in >> n;
	if (!in.good())
	  return err("Expected integer after 'number-of-variables:' in "
		     "in HEADER section");
      } else if (tk == "number-of-constraints:") {
	in >> m;
	if (!in.good())
	  return err("Expected integer after 'number-of-variables:' in "
		     "in HEADER section");
      } else if (tk == "default-bound:") {
	if (nt == YET_UNKOWN_NUMBER_TYPE)
	  return err("'default-bound:' must be preceded by 'number-type:' "
		     "in HEADER section");
	std::string bnd = token();
	if (bnd == "nonnegative")
	  ; // no need to do anything at all as this is the default
	else if (bnd == "free")
	  default_fl = false;
	else
	  return err1("Unknown bound type '%' in 'default-bound:' clause "
		      "in HEADER section", bnd);
      } else 
	return err1("Unknown tag '%' in HEADER section", tk);
      
      tk = token();
    }

    if (nt == YET_UNKOWN_NUMBER_TYPE)
      return err("No 'number-type:' given in HEADER");
    if (n < 0)
      return err("No 'number-of-variables:' given in HEADER");
    if (m < 0)
      return err("No 'number-of-constraints:' given in HEADER");
    return true;
  }

  // Reads the constraints, returning false in case of an error.
  // Note: this method sets error_message_ appropriately in case of errors.
  bool read_CONSTRAINTS() {
    // Note: "CONSTRAINTS" has already been parsed at this point.
    
    ANT val;
    for (int i = 0; i < m; ++i) {
      // read row of A:
      for (int j = 0; j < n; ++j) {
	if (!read_number(val))
	  return err2("Expected number in CONSTRAINTS section at A[%,%]",
		      tostr(i), tostr(j));
	A_[j][i] = val;
      }
    
      // read row-type:
      std::string rt = token();
      if (rt == "<=")
	row_type_[i] = LESS_EQUAL;
      else if (rt == "=" || rt == "==")
	row_type_[i] = EQUAL;
      else if (rt == ">=")
	row_type_[i] = GREATER_EQUAL;
      else
	return err2("Expected one of {'<=','=','>='} but found '%' in "
		    "row % of CONSTRAINTS section", rt, tostr(i));

      // read b-entry:
      if (!read_number(val))
	return err1("Expected number in CONSTRAINTS section at b[%]", 
		    tostr(i));
      b_[i] = val;
    }

    return true;
  }

  // Reads the D-matrix, if need be, returning false in case of an error.
  // Note: this method sets error_message_ appropriately in case of errors.
  // Note: this method reads the next token into tk.
  bool read_D_MATRIX(std::string& tk) {
    // check if there is a D-MATRIX section at all:
    tk = token();
    if (tk != "D-MATRIX")
      return true;

    ANT val;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
	if (!read_number(val))
	  return err2("Expected number in D-MATRIX section at D[%,%]",
		      tostr(i), tostr(j));
	D_[i][j] = val;
      }
    }

    tk = token();
    return true;
  }

  // Reads the c-vector and the bounds, returning false in case of an error.
  // Note: this method sets error_message_ appropriately in case of errors.
  bool read_C_VECTOR_AND_BOUNDS(const std::string& tk) {
    // safety check:
    if (tk != "C-VECTOR-AND-BOUNDS")
      return err1("Expected C-VECTOR-AND-BOUNDS but found '%'", tk);

    ANT val;
    for (int i = 0; i < n; ++i) {
      // read c-entry:
      if (!read_number(val))
	return err1("Expected number in C-VECTOR-AND-BOUNDS section for "
		    "variable %", tostr(i));
      c_[i] = val;
      
      // read bounds, if need be:
      whitespace();
      if (in.peek() == '(') {
	in.get();
	whitespace();
	if (in.peek() == ',') {
	  in.get();
	  fl_[i] = false;
	} else {
	  if (!read_number(val))
	    return err1("Expected number in C-VECTOR-AND-BOUNDS section for "
			"finite lower bound of variable %", tostr(i));
	  l_[i] = val;
	  whitespace();
	  if (in.get() != ',')
	    return err1("Expected ',' in C-VECTOR-AND-BOUNDS section after "
			"finite lower bound of variable %", tostr(i));
	}
	whitespace();
	if (in.peek() == ')') {
	  in.get();
	  fu_[i] = false;
	} else {
	  if (!read_number(val))
	    return err1("Expected number in C-VECTOR-AND-BOUNDS section for "
			"finite upper bound of variable %", tostr(i));
	  u_[i] = val;
	  whitespace();
	  if (in.get() != ')')
	    return err1("Expected ')' in C-VECTOR-AND-BOUNDS section after "
			"finite upper bound of variable %", tostr(i));
	}
      }
    }

    return true;
  }
};

int main(int argnr, char **argv)
{
  // output usage information:
  if (argnr > 1) {
    std::cerr << "Usage: " << argv[0] << " < datafile.data\n\n"
	      << "The resulting MPS-file is written to standard-out, errors "
	      << "go to standard-error.\n\n"
	      << "See data2mps.cin for an example 'datafile.data' file.\n";
    return 1;
  }

  // read data-file from standard in:
  Data_reader data(std::cin);
  if (!data.good()) {
    std::cerr << "An error occurred while parsing the input data file: \n"
	      << data.error_message() << "\n";
    return 2;
  }

  // write MPS-file to standard out:

}
