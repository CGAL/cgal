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
// file          : test/_QP_solver/test_QP_solver.C
// package       : $CGAL_Package: _QP_solver $
//
// revision      : 0.2
// revision_date : 2000/08/21
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
//                 Frans Wessendorp <fransw@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test program for the QP solver
// ============================================================================

#include <CGAL/QP_solver/gmp_double.h>
#include <CGAL/QP_solver.h>
#include <CGAL/QP_full_exact_pricing.h>
#include <CGAL/QP_partial_exact_pricing.h>
#include <CGAL/QP_full_filtered_pricing.h>
#include <CGAL/QP_partial_filtered_pricing.h>
#include <CGAL/QP_solver/Double.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#include <cassert>
#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <locale>

// We first define two template classes, Vector<ET> and Matrix<ET>, to 
// hold the input matrices and vector.
template<typename ET_>
class Vector : public std::vector<ET_>
{
  typedef std::vector<ET_> base;
  
public:
  Vector(typename base::size_type i) : base(i) {}
  Vector() : base() {}
};

template<typename ET_>
class Matrix : public std::vector< Vector<ET_> >
{
};

// Next, we need a functor that
template<typename T_>
struct Begin
  : public std::unary_function< Vector<T_>,
				typename Vector<T_>::const_iterator > {
  typedef typename Vector<T_>::const_iterator result_type;
  result_type operator () ( const Vector<T_>& v) const { return v.begin(); }
};

template<typename ET_,
	 typename IT_,
	 typename Is_linear_,
	 typename Is_symmetric_,
	 typename Has_equalities_only_and_full_rank_,
	 typename Is_in_standard_form_>
struct Rep {
  typedef typename CGAL::Join_input_iterator_1< 
    typename Matrix<IT_>::const_iterator, Begin<IT_> > Vector_iterator;
  typedef typename Vector<IT_>::const_iterator Entry_iterator;
  typedef IT_ Input_type;
  
  
  typedef  Vector_iterator  A_iterator;
  typedef   Entry_iterator  B_iterator;
  typedef   Entry_iterator  C_iterator;
  typedef  Vector_iterator  D_iterator;
  typedef   Entry_iterator  L_iterator;
  typedef   Entry_iterator  U_iterator;
  typedef   std::vector<bool>::const_iterator  FL_iterator;
  typedef   std::vector<bool>::const_iterator  FU_iterator;

  enum Row_type { LESS_EQUAL = -1, EQUAL, GREATER_EQUAL};
  typedef  Row_type*  Row_type_iterator;
  //    typedef  CGAL::QP_const_value_iterator<Row_type>  Row_type_iterator;
    
  //typedef  CGAL::Double  ET;
  typedef ET_  ET;
  typedef Is_linear_  Is_linear;
  typedef Is_symmetric_  Is_symmetric;
  typedef Has_equalities_only_and_full_rank_  Has_equalities_only_and_full_rank;
  typedef Is_in_standard_form_ Is_in_standard_form;
};

template<typename Rep>
void init_row_types(std::vector<int>& rel,typename Rep::Row_type*& row_types);

template<typename T>
bool read_explicit_bounds(std::ifstream& from, std::vector<bool>& fl,
            Vector<T>& l, std::vector<bool>& fu, Vector<T>& u, int n,
            CGAL::Tag_true Is_in_standard_form);

template<typename T>
bool read_explicit_bounds(std::ifstream& from, std::vector<bool>& fl,
            Vector<T>& l, std::vector<bool>& fu, Vector<T>& u, int n,
            CGAL::Tag_false Is_in_standard_form);

template<typename T>
bool read_instance(std::ifstream& from, Matrix<T>& a,
		   std::vector<int>& rel, Vector<T>& b,
		   Vector<T>& c, Matrix<T>& d, int& n, CGAL::Tag_true Is_linear);

template<typename T>
bool read_instance(std::ifstream& from, Matrix<T>& a,
		   std::vector<int>& rel, Vector<T>& b,
		   Vector<T>& c, Matrix<T>& d, int& n, CGAL::Tag_false Is_linear);
	
template<typename T>
bool read_matrix(std::ifstream& from, int m, int n,
		 //typename Vector<T>::size_type m,
		 //typename Vector<T>::size_type n,
		 Matrix<T>& mat);
	
	
template<typename T>
bool read_matrix_transposed(std::ifstream& from, int m, int n,
			    //typename Vector<T>::size_type m,
			    //typename Vector<T>::size_type n,
			    Matrix<T>& mat);
	

template<typename T>
bool read_vector(std::ifstream& from, int length,
		 Vector<T>& vec);

bool read_num_entry(std::ifstream& from, CGAL::Gmpq& entry);

bool read_num_entry(std::ifstream& from, CGAL::Double& entry);

bool read_num_entry(std::ifstream& from, CGAL::Gmpz& entry);

bool read_int_vector(std::ifstream& from, int length, std::vector<int>& vec);

bool read_bool_vector(std::ifstream& from, int length, std::vector<bool>& vec);

void read_ws(std::istream& from);

bool read_token(std::ifstream& from, const std::string& block_del);
	
void init_tag_names_table();
	
bool read_tags(std::ifstream& from, 
	       std::map<std::string, char>& read_names_table);	

bool map_tags(std::ifstream& from, int verbose, int pricing_strategy_index,
	      const std::map<std::string, char>& read_names_table);
	
bool insert_name(const std::string& name, const char& value,
		 std::map<std::string, char>& read_names_table);
	
void error(const std::string& msg);

void message(const std::string& msg, bool mode);

void init_pricing_strat_abrev_table();

int map_pricing_strategy_abrev(const std::string& abbrev);

template <typename Rep>
void set_pricing_strategy(CGAL::QP_pricing_strategy<Rep>*& strat, int index);


std::vector<std::string> tag_names_table(0);
std::vector<std::string> pricing_strat_abrev_table(0);

const bool debug = false;

template<typename ET>
struct Input_type_selector {};

template<>
struct Input_type_selector<CGAL::Gmpq> {
  typedef CGAL::Gmpq Input_type;
};

template<>
struct Input_type_selector<CGAL::Double> {
  typedef double Input_type;
};

template<>
struct Input_type_selector<CGAL::Gmpz> {
  typedef CGAL::Gmpz Input_type;
};

template<typename ET,
	 typename Is_in_standard_form,
	 typename Is_linear,
	 typename Is_symmetric,
	 typename Has_equalities_only_and_full_rank>
bool doIt(int verbose, int pricing_strategy_index, std::ifstream& from);  
 

int main( int argc, char** argv) {
  init_tag_names_table();
  init_pricing_strat_abrev_table();
  int verbose = 0;
  unsigned int pricing_strategy_index = 0; // to prevent GCC from
  // outputing a warning
  std::map<std::string, char> read_names_table;

  const bool readFromStdIn = argc <= 3;
  if (readFromStdIn)
    std::cout << "reading from standard in..." << std::endl;
  else {
    verbose = atoi(argv[1]);
    pricing_strategy_index = map_pricing_strategy_abrev(argv[2]);
    if (pricing_strategy_index == pricing_strat_abrev_table.size()) {
      std::cout << "unidentified pricing strategy" << std::endl;
      return 1;
    }
  }
	
  bool success = true;
  for(int i=3; i<argc || (readFromStdIn && !std::cin.eof()); ++i) {

    // read verbosity and pricing-strategy
    if (readFromStdIn) {
      std::cin >> verbose;
      std::string strategy;
      std::cin >> strategy;
      pricing_strategy_index = map_pricing_strategy_abrev(strategy);
      if (pricing_strategy_index == pricing_strat_abrev_table.size()) {
	std::cout << "unidentified pricing strategy" << std::endl;
	return 1;
      }
    }
	  
    std::string fileName;
    if (readFromStdIn)
      std::cin >> fileName;
    else
      fileName = std::string(argv[i]);
    std::ifstream from(fileName.c_str());
    if (!from) {
      std::cout << "could not open file: " << fileName << "\n";
    } else {
      std::cout << "processing file: " << fileName << "\n";
      read_tags(from, read_names_table);
      if (!map_tags(from, verbose, pricing_strategy_index,
		    read_names_table))
	success = false;
      from.close();
    }
    if (readFromStdIn)
      read_ws(std::cin);
  }

  if (!success)
    std::cout << "Warning: some test cases were not handled successfully." << std::endl;
  return success? 0 : 2;
}


template<typename ET,
	 typename Is_in_standard_form,
	 typename Is_linear,
	 typename Is_symmetric,
	 typename Has_equalities_only_and_full_rank>
bool doIt(int verbose, int pricing_strategy_index, std::ifstream& from) {
  typedef typename Input_type_selector<ET>::Input_type Input_type;
  typedef Rep<ET, Input_type, Is_linear, Is_symmetric,
    Has_equalities_only_and_full_rank, Is_in_standard_form> Repr;
    
  Matrix<Input_type>  A;
  typename Repr::Row_type* row_types;
  Vector<Input_type>  b;
  Vector<Input_type>  c;
  Matrix<Input_type>  D;
  std::vector<int> rel;
  std::vector<bool> fl, fu;
  Vector<Input_type> l,u;
  int n;
  
  CGAL::QP_pricing_strategy<Repr>*
    strat(static_cast<CGAL::QP_pricing_strategy<Repr>*>(0));
  bool instance_read, sol_solver_valid;	


  instance_read = read_instance<Input_type>(from, A, rel, b, c, D, n,
					    Is_linear());
				
  instance_read = instance_read && read_explicit_bounds(from, fl, l, fu, u, n,
                        Is_in_standard_form());

  if (instance_read) {
    init_row_types<Repr>(rel, row_types);
    set_pricing_strategy<Repr>(strat, pricing_strategy_index);	
    CGAL::QP_solver< Repr > solver( A.size(), rel.size(),
		typename Repr::Vector_iterator( A.begin()),
		b.begin(), c.begin(),
		typename Repr::Vector_iterator( D.begin()),
		row_types, fl.begin(), l.begin(), fu.begin(), u.begin(), *strat);
    solver.set_verbosity( verbose);
    sol_solver_valid = solver.is_solution_valid();
    std::cout << "valid: " << sol_solver_valid << std::endl;
    delete strat;
    delete[] row_types;
    return sol_solver_valid;
  } else {
    error("could not read problem instance");
    return false;
  }
}
 


template<typename Rep>
void init_row_types(std::vector<int>& rel,typename Rep::Row_type*& row_types) {
  int size = rel.size();
  //delete[] row_types;
  row_types = new typename Rep::Row_type[size];
	
  for (int i = 0; i < size; ++i) {
    row_types[i] = typename Rep::Row_type(rel[i]); 		
  }
	
}

template<typename T>
bool read_explicit_bounds(std::ifstream& from, std::vector<bool>& fl,
            Vector<T>& l, std::vector<bool>& fu, Vector<T>& u, int n,
            CGAL::Tag_true Is_in_standard_form) {
  return true;
}

template<typename T>
bool read_explicit_bounds(std::ifstream& from, std::vector<bool>& fl,
            Vector<T>& l, std::vector<bool>& fu, Vector<T>& u, int n,
            CGAL::Tag_false Is_in_standard_form) {

  bool token_read, matrix_read;
  
  // read lower bound finiteness vector fl
  token_read = read_token(from, "fl");
  if (token_read) {
    message("fl token matched", debug);
    matrix_read = read_bool_vector(from, n, fl);
    if (!matrix_read) {
      error("could not read matrix fl");
      return false;
    }

  } else {
    error("could not match token fl");
    return false;
  }

  // read finite lower bounds l
  token_read = read_token(from, "l");
  if (token_read) {
    message("l token matched", debug);
    matrix_read = read_vector<T>(from, n, l);
    if (!matrix_read) {
      error("could not read matrix l");
      return false;
    }

  } else {
    error("could not match token l");
    return false;
  }

  // read upper bound finiteness vector fu 
  token_read = read_token(from, "fu");
  if (token_read) {
    message("fu token matched", debug);
    matrix_read = read_bool_vector(from, n, fu);
    if (!matrix_read) {
      error("could not read matrix fu");
      return false;
    }

  } else {
    error("could not match token fu");
    return false;
  }

  // read finite upper bounds u
  token_read = read_token(from, "u");
  if (token_read) {
    message("u token matched", debug);
    matrix_read = read_vector<T>(from, n, u);
    if (!matrix_read) {
      error("could not read matrix u");
      return false;
    }

  } else {
    error("could not match token u");
    return false;
  }
  return true;
}
          
    


template<typename T>
bool read_instance(std::ifstream& from, Matrix<T>& a,
		   std::vector<int>& rel, Vector<T>& b,
		   Vector<T>& c, Matrix<T>& d, int& n, CGAL::Tag_false) {	//QP case
	
  int m;
  bool token_read, matrix_read;
		
  // read dimensions
  token_read = read_token(from, "dimensions");
  if (token_read) {
    message("dimensions token matched", debug);
    from >> m; from >> n;
    if (!from.good()) {
      error("could not read dimensions");
      return false;
    }
    std::cout << "m= " << m << " n= " << n << "\n";
  } else {
    error("could not match token dimensions");
    return false;
  }
	
  // read constraint matrix A
  token_read = read_token(from, "A");
  if (token_read) {
    message("A token matched", debug);
    matrix_read = read_matrix_transposed<T>(from, m, n, a);
    if (!matrix_read) {
      error("could not read matrix A");
      return false;
    }
  } else {
    error("could not match token A");
    return false;
  }
	
  // read rel operators
  token_read = read_token(from, "rel");
  if (token_read) {
    message("rel token matched", debug);
    matrix_read = read_int_vector(from, m, rel);
    if (!matrix_read) {
      error("could not read matrix rel");
      return false;
    }
  } else {
    error("could not match token rel");
    return false;
  }
	
  // read rhs b
  token_read = read_token(from, "b");
  if (token_read) {
    message("b token matched", debug);
    matrix_read = read_vector<T>(from, m, b);
    if (!matrix_read) {
      error("could not read matrix b");
      return false;
    }

  } else {
    error("could not match token b");
    return false;
  }
	
  // read objfunc c
  token_read = read_token(from, "c");
  if (token_read) {
    message("c token matched", debug);
    matrix_read = read_vector<T>(from, n, c);
    if (!matrix_read) {
      error("could not read matrix c");
      return false;
    }

  } else {
    error("could not match token c");
    return false;
  }
  // read objfunc D
  token_read = read_token(from, "D");
  if (token_read) {
    message("D token matched", debug);
    matrix_read = read_matrix<T>(from, n, n, d);
    if (!matrix_read) {
      error("could not read matrix D");
      return false;
    }
  } else {
    error("could not match token D");
    return false;
  }
  return true;
}

template<typename T>
bool read_instance(std::ifstream& from, Matrix<T>& a,
		   std::vector<int>& rel, Vector<T>& b,
		   Vector<T>& c, Matrix<T>& d, int& n, CGAL::Tag_true) { //LP case
	
  int m;
  bool token_read, matrix_read;
		
  // read dimensions
  token_read = read_token(from, "dimensions");
  if (token_read) {
    message("dimensions token matched", debug);
    from >> m; from >> n;
    if (!from.good()) {
      error("could not read dimensions");
      return false;
    }
    std::cout << "m= " << m << " n= " << n << "\n";
  } else {
    error("could not match token dimensions");
    return false;
  }
	
  // read constraint matrix A
  token_read = read_token(from, "A");
  if (token_read) {
    message("A token matched", debug);
    matrix_read = read_matrix_transposed<T>(from, m, n, a);
    if (!matrix_read) {
      error("could not read matrix A");
      return false;
    }
  } else {
    error("could not match token A");
    return false;
  }
	
  // read rel operators
  token_read = read_token(from, "rel");
  if (token_read) {
    message("rel token matched", debug);
    matrix_read = read_int_vector(from, m, rel);
    if (!matrix_read) {
      error("could not read matrix rel");
      return false;
    }
  } else {
    error("could not match token rel");
    return false;
  }
	
  // read rhs b
  token_read = read_token(from, "b");
  if (token_read) {
    message("b token matched", debug);
    matrix_read = read_vector<T>(from, m, b);
    if (!matrix_read) {
      error("could not read matrix b");
      return false;
    }

  } else {
    error("could not match token b");
    return false;
  }
	
  // read objfunc c
  token_read = read_token(from, "c");
  if (token_read) {
    message("c token matched", debug);
    matrix_read = read_vector<T>(from, n, c);
    if (!matrix_read) {
      error("could not read matrix c");
      return false;
    }

  } else {
    error("could not match token c");
    return false;
  }
  return true;
}


template<typename T>
bool read_matrix(std::ifstream& from, int m, int n,
		 //typename Vector<T>::size_type m,
		 //typename Vector<T>::size_type n,
		 Matrix<T>& mat) {
	
  std::ostringstream msgstr; 
  bool entry_read = true;	
  mat.assign( m, Vector<T>(n));
  for (int row = 0; row < m; ++row) {
    for (int col = 0; col < n; ++col) {
      entry_read = read_num_entry(from, mat[row][col]);
      if (!entry_read) {
	return false;
      }
      msgstr.str("");
      msgstr << "row: " << row << " col: " << col
	     << " entry: " << mat[row][col];
      message(msgstr.str(), debug);
    }
  } 
  return true;
}

template<typename T>
bool read_matrix_transposed(std::ifstream& from, int m, int n,
			    //typename Vector<T>::size_type m,
			    //typename Vector<T>::size_type n,
			    Matrix<T>& mat) {
	
  std::ostringstream msgstr; 
  bool entry_read = true;
  mat.assign( n, Vector<T>(m));
  for (int row = 0; row < m; ++row) {
    for (int col = 0; col < n; ++col) {
      entry_read = read_num_entry(from, mat[col][row]);
      if (!entry_read) {
	return false;
      }
      msgstr.str("");
      msgstr << "row: " << row << " col: " 
	     << col << " entry: " << mat[col][row];
      message(msgstr.str(), debug);
    }
  } 
  return true;
}

template<typename T>
bool read_vector(std::ifstream& from, int length, Vector<T>& vec) {
	
  std::ostringstream msgstr;
  bool entry_read = true;
  vec.assign(length, T());
  for (int i = 0; i < length; ++i) {
    entry_read = read_num_entry(from, vec[i]);
    if (!entry_read) {
      return false;
    }
    msgstr.str("");
    msgstr << "row/col: " << i << " entry: " << vec[i];
    message(msgstr.str(), debug);
  }
  return true;
}

bool read_num_entry(std::ifstream& from, CGAL::Gmpq& entry) {
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

bool read_num_entry(std::ifstream& from, double& entry) {
  from >> entry;
  return from.good();
}

bool read_num_entry(std::ifstream& from, CGAL::Gmpz& entry) {
  from >> entry;
  return from.good();
}

	
bool read_int_vector(std::ifstream& from, int length, std::vector<int>& vec) {
	
  std::ostringstream msgstr;
  int i = 0;
  vec.assign(length, 0);
  while ((i < length) && from.good()) {
    from >> vec[i];
    msgstr.str("");
    msgstr << "row/col: " << i << " entry: " << vec[i];
    message(msgstr.str(), debug);
    i = i + 1;
  }
  return from.good();
}

bool read_bool_vector(std::ifstream& from, int length, std::vector<bool>& vec) {
  std::ostringstream msgstr;
   bool b;
  int i = 0;
  vec.assign(length, 0);
  while ((i < length) && from.good()) {
    // the following two statements are a workaround due to implicit C++
    // optimization crap
    from >> b;
    vec[i] = b;
    msgstr.str("");
    msgstr << "row/col: " << i << " entry: " << vec[i];
    message(msgstr.str(), debug);
    i = i + 1;
  }
  return from.good();
} 

bool read_tags(std::ifstream& from,
	       std::map<std::string, char>& read_names_table){
	
  std::ostringstream msgstr;
  std::string tag_name, dummy;
  char ch, tag_value;
  bool ok, unique;
  unsigned int i = 0;
	
  //read comments' section
  read_ws(from);
  from.get(ch);
  while (ch == '#') {
    std::getline(from, dummy);
    read_ws(from);
    from.get(ch);
  }
  from.putback(ch);
	
  //read leading whitespace
  ok = read_token(from, "tags");
  if (ok) {
    message("tags token matched", debug);
  } else {
    error("could not match tags token");
    return false;
  }
  unique = true;
  read_names_table.clear();
  while ( (i < tag_names_table.size()) && ok && unique && from.good()) {
    from >> tag_name;
    from >> tag_value;
    msgstr.str("");
    msgstr << "entry " << i << " : " << std::endl;
    msgstr << "tag name: " << tag_name << std::endl;
    msgstr << "tag value: " << tag_value;
    message(msgstr.str(), debug);
    unique = insert_name(tag_name, tag_value, read_names_table);
    i = i + 1;
  }
  if ((i == tag_names_table.size()) && ok && unique && from.good()) {	
    return true;
  } else {
    error("error in read tags");
    return false;
  }
}

bool map_tags(std::ifstream& from, int verbose, int pricing_strategy_index,
	      const std::map<std::string, char>& read_names_table) {
  int offset;
  std::map<std::string, char>::const_iterator p;
  int index = 0;
  offset = tag_names_table.size() - 2; // we handle is_in_standard_form 
                                       // seperately
  for (unsigned int i = 0; i < tag_names_table.size(); ++i) {
    p = read_names_table.find(tag_names_table[i]);
    if (p != read_names_table.end()) {
      if (i==0) { // r,f,i
	switch (p->second) {
	case 'r':	index = index + 0;
	  break;
	case 'f':	index = index + 8;
	  break;
	case 'i':	index = index + 16;
	  break;
	default : 
	  error("impossible tag value");
	  exit(-1);
	}
      } else if (i==4) {  // is_in_standard_form
	if (p->second == '1')
	  index += 24;
      } else { // 0,1
	switch (p->second) {
	case '0':	;
	  break;
	case '1':	index += (1 << (offset - i));
	  break;
	default :
	  error("impossible tag value");
	  exit(-1);
	}
      }
    } else {
      error("unspecified tags");
      exit(-1);
    }
  }

  using CGAL::Tag_false;
  using CGAL::Tag_true;

  bool success;
  switch (index) {
  case  0: 	
    success = doIt<CGAL::Gmpq,Tag_false,Tag_false,Tag_false,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case  1:
    success = doIt<CGAL::Gmpq,Tag_false,Tag_false,Tag_false,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    
  case  2:	
    success = doIt<CGAL::Gmpq,Tag_false,Tag_false,Tag_true,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case  3:	
    success = doIt<CGAL::Gmpq,Tag_false,Tag_false,Tag_true,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    
    // next 2 cases can be excluded by assuming is_Symmetric == 1 for
    // LPs
    //case  4:	
    // success = doIt<CGAL::Gmpq,Tag_false,Tag_true,Tag_false,
    //			Tag_false>
    //			(verbose, pricing_strategy_index, from);
    //		break;
    //case  5:	
    //  success = doIt<CGAL::Gmpq,Tag_false,Tag_true,Tag_false,
    //			Tag_true>
    //			(verbose, pricing_strategy_index, from);
    //		break;
  case  6:	
    success = doIt<CGAL::Gmpq,Tag_false,Tag_true,Tag_true,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case  7:	
    success = doIt<CGAL::Gmpq,Tag_false,Tag_true,Tag_true,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
  case  8:	
    success = doIt<CGAL::Double,Tag_false,Tag_false,Tag_false,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case  9:	
    success = doIt<CGAL::Double,Tag_false,Tag_false,Tag_false,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
  case 10:	
    success = doIt<CGAL::Double,Tag_false,Tag_false,Tag_true,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 11:	
    success = doIt<CGAL::Double,Tag_false,Tag_false,Tag_true,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    // next 2 cases can be excluded by assuming is_Symmetric == 1 for
    // LPs	
    //case 12:	
    //  success = doIt<CGAL::Double,Tag_false,Tag_true,Tag_false,
    //			Tag_false>
    //			(verbose, pricing_strategy_index, from);
    //		break;
    //case 13:	
    //  success = doIt<CGAL::Double,Tag_false,Tag_true,Tag_false,
    //			Tag_true>
    //			(verbose, pricing_strategy_index, from);
    //		break;
  case 14:	
    success = doIt<CGAL::Double,Tag_false,Tag_true,Tag_true,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 15:	
    success = doIt<CGAL::Double,Tag_false,Tag_true,Tag_true,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 16:	
    success = doIt<CGAL::Gmpz,Tag_false,Tag_false,Tag_false,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 17:	
    success = doIt<CGAL::Gmpz,Tag_false,Tag_false,Tag_false,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 18:	
    success = doIt<CGAL::Gmpz,Tag_false,Tag_false,Tag_true,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 19:	
    success = doIt<CGAL::Gmpz,Tag_false,Tag_false,Tag_true,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    // next 2 cases can be excluded by assuming is_Symmetric == 1 for
    // LPs
    //case 20:	
    // success = doIt<CGAL::Gmpz,Tag_false,Tag_true,Tag_false,
    //			Tag_false>
    //			(verbose, pricing_strategy_index, from);
    //		break;
    //case 21:	
    //  success = doIt<CGAL::Gmpz,Tag_false,Tag_true,Tag_false,
    //			Tag_true>
    //			(verbose, pricing_strategy_index, from);
    //		break;
  case 22:	
    success = doIt<CGAL::Gmpz,Tag_false,Tag_true,Tag_true,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 23:	
    success = doIt<CGAL::Gmpz,Tag_false,Tag_true,Tag_true,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    
  case  0+24: 	
    success = doIt<CGAL::Gmpq,Tag_true,Tag_false,Tag_false,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case  1+24:
    success = doIt<CGAL::Gmpq,Tag_true,Tag_false,Tag_false,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    
  case  2+24:	
    success = doIt<CGAL::Gmpq,Tag_true,Tag_false,Tag_true,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case  3+24:	
    success = doIt<CGAL::Gmpq,Tag_true,Tag_false,Tag_true,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    
    // next 2 cases can be excluded by assuming is_Symmetric == 1 for
    // LPs
    //case  4+24:	
    // success = doIt<CGAL::Gmpq,Tag_true,Tag_true,Tag_false,
    //			Tag_false>
    //			(verbose, pricing_strategy_index, from);
    //		break;
    //case  5+24:	
    //  success = doIt<CGAL::Gmpq,Tag_true,Tag_true,Tag_false,
    //			Tag_true>
    //			(verbose, pricing_strategy_index, from);
    //		break;
  case  6+24:	
    success = doIt<CGAL::Gmpq,Tag_true,Tag_true,Tag_true,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case  7+24:	
    success = doIt<CGAL::Gmpq,Tag_true,Tag_true,Tag_true,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
  case  8+24:	
    success = doIt<CGAL::Double,Tag_true,Tag_false,Tag_false,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case  9+24:	
    success = doIt<CGAL::Double,Tag_true,Tag_false,Tag_false,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
  case 10+24:	
    success = doIt<CGAL::Double,Tag_true,Tag_false,Tag_true,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 11+24:	
    success = doIt<CGAL::Double,Tag_true,Tag_false,Tag_true,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    // next 2 cases can be excluded by assuming is_Symmetric == 1 for
    // LPs	
    //case 12+24:	
    //  success = doIt<CGAL::Double,Tag_true,Tag_true,Tag_false,
    //			Tag_false>
    //			(verbose, pricing_strategy_index, from);
    //		break;
    //case 13+24:	
    //  success = doIt<CGAL::Double,Tag_true,Tag_true,Tag_false,
    //			Tag_true>
    //			(verbose, pricing_strategy_index, from);
    //		break;
  case 14+24:	
    success = doIt<CGAL::Double,Tag_true,Tag_true,Tag_true,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 15+24:	
    success = doIt<CGAL::Double,Tag_true,Tag_true,Tag_true,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 16+24:	
    success = doIt<CGAL::Gmpz,Tag_true,Tag_false,Tag_false,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 17+24:	
    success = doIt<CGAL::Gmpz,Tag_true,Tag_false,Tag_false,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 18+24:	
    success = doIt<CGAL::Gmpz,Tag_true,Tag_false,Tag_true,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 19+24:	
    success = doIt<CGAL::Gmpz,Tag_true,Tag_false,Tag_true,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    // next 2 cases can be excluded by assuming is_Symmetric == 1 for
    // LPs
    //case 20+24:	
    // success = doIt<CGAL::Gmpz,Tag_true,Tag_true,Tag_false,
    //			Tag_false>
    //			(verbose, pricing_strategy_index, from);
    //		break;
    //case 21+24:	
    //  success = doIt<CGAL::Gmpz,Tag_true,Tag_true,Tag_false,
    //			Tag_true>
    //			(verbose, pricing_strategy_index, from);
    //		break;
  case 22+24:	
    success = doIt<CGAL::Gmpz,Tag_true,Tag_true,Tag_true,
      Tag_false>
      (verbose, pricing_strategy_index, from);
    break;
    
  case 23+24:	
    success = doIt<CGAL::Gmpz,Tag_true,Tag_true,Tag_true,
      Tag_true>
      (verbose, pricing_strategy_index, from);
    break;
    
  default:	
    std::cout << "map_tags: unidentified case" << std::endl;
    success = false;
  }

  return success;
}

void init_tag_names_table() {
  tag_names_table.clear();
  tag_names_table.push_back("input_data_type");
  tag_names_table.push_back("is_linear");
  tag_names_table.push_back("is_symmetric");
  tag_names_table.push_back("has_equalities_only_and_full_rank");	
  tag_names_table.push_back("is_in_standard_form");	
}


bool insert_name(const std::string& name, const char& value,
		 std::map<std::string, char>& read_names_table) {
	
  std::pair<std::string, char> tmp(name, value);
  std::pair<std::map<std::string, char>::iterator,
    bool> p = read_names_table.insert(tmp);
  return p.second;
}

 
void read_ws(std::istream& from) {
  char c;
  while (from.get(c)) {
    if (!isspace(c)) {
      from.putback(c);
      break;
    }
  }
}	 


bool read_token(std::ifstream& from, const std::string& block_del) {
	
  std::string token;
  from >> token;
  if (token == block_del) {
    return true;
  } else {
    return false;
  }
}

void error(const std::string& msg) {
  std::cerr << msg << std::endl;
}

void message(const std::string& msg, bool mode) {
  if (mode) {
    std::cout << msg << std::endl;
  }
}

void init_pricing_strat_abrev_table() {
  pricing_strat_abrev_table.clear();
  pricing_strat_abrev_table.push_back("fe");
  pricing_strat_abrev_table.push_back("ff");
  pricing_strat_abrev_table.push_back("pe");
  pricing_strat_abrev_table.push_back("pf");
}	

int map_pricing_strategy_abrev(const std::string& abrev) {
  int ind = pricing_strat_abrev_table.size();
  for (unsigned int i = 0; i < pricing_strat_abrev_table.size(); ++i) {
    if (abrev == pricing_strat_abrev_table[i]) {
      ind = i;
    }
  }
  return ind;
}
template <typename Rep>
void set_pricing_strategy(CGAL::QP_pricing_strategy<Rep>*& strat, int index) {
  switch (index) {
  case 0	:	strat = new CGAL::QP_full_exact_pricing<Rep>;
    break;
  case 1	:	strat = new CGAL::QP_full_filtered_pricing<Rep,
					    typename Rep::Input_type>;
    break;
  case 2	:	strat = new CGAL::QP_partial_exact_pricing<Rep>;
    break;
  case 3	:	strat = new CGAL::QP_partial_filtered_pricing<Rep,
					    typename Rep::Input_type>;
    break;
  }
}


// ===== EOF ==================================================================
