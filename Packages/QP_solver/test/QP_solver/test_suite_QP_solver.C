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
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test program for the QP solver
// ============================================================================

#include <CGAL/QPE_solver.h>
#include <CGAL/QPE_full_exact_pricing.h>
#include <CGAL/QPE_partial_exact_pricing.h>
#include <CGAL/QPE_full_filtered_pricing.h>
#include <CGAL/QPE_partial_filtered_pricing.h>
#include <CGAL/_QP_solver/Double.h>
//#include <CGAL/Quotient.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#include <cassert>
#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <locale>

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


template<typename T_>
struct Begin
    : public CGAL_STD::unary_function< Vector<T_>, typename Vector<T_>::const_iterator > {
    typedef typename Vector<T_>::const_iterator result_type;
	result_type operator () ( const Vector<T_>& v) const { return v.begin(); }
};

template<typename ET_,
	typename IT_,
	typename Is_linear_,
 	typename Is_symmetric_,
	typename Has_no_inequalities_>
struct Rep {
	typedef typename CGAL::QPE_transform_iterator_1< typename Matrix<IT_>::const_iterator, Begin<IT_> > Vector_iterator;
	typedef typename Vector<IT_>::const_iterator Entry_iterator;


    typedef  Vector_iterator  A_iterator;
    typedef   Entry_iterator  B_iterator;
    typedef   Entry_iterator  C_iterator;
    typedef  Vector_iterator  D_iterator;

    enum Row_type { LESS_EQUAL = -1, EQUAL, GREATER_EQUAL};
    typedef  Row_type*  Row_type_iterator;
    //    typedef  CGAL::QPE_const_value_iterator<Row_type>  Row_type_iterator;
    
    //typedef  GMP::Double  ET;
    typedef ET_  ET;
    typedef  Is_linear_  Is_linear;
    typedef  Is_symmetric_  Is_symmetric;
    typedef  Has_no_inequalities_  Has_no_inequalities;
};


//struct Rep {
    //typedef  GMP::Double  ET;
//    typedef CGAL::Gmpq  ET;
//    typedef  Vector_iterator  A_iterator;
//    typedef   Entry_iterator  B_iterator;
//    typedef   Entry_iterator  C_iterator;
//    typedef  Vector_iterator  D_iterator;

//    enum Row_type { LESS_EQUAL = -1, EQUAL, GREATER_EQUAL};
//    typedef  Row_type*  Row_type_iterator;
    //    typedef  CGAL::QPE_const_value_iterator<Row_type>  Row_type_iterator;
    

//    typedef  CGAL::Tag_false  Is_linear;
//    typedef  CGAL::Tag_true  Is_symmetric;
//    typedef  CGAL::Tag_true  Has_no_inequalities;
//};

//typedef std::map<std::string, char> tag_names_values_map;

template<typename Rep>
void init_row_types(std::vector<int>& rel,typename Rep::Row_type*& row_types);

template<typename T>
bool read_instance(std::ifstream& from, Matrix<T>& a,
	std::vector<int>& rel, Vector<T>& b,
	Vector<T>& c, Matrix<T>& d);

	
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

bool read_num_entry(std::ifstream& from, GMP::Double& entry);

bool read_num_entry(std::ifstream& from, CGAL::Gmpz& entry);

bool read_int_vector(std::ifstream& from, int length, std::vector<int>& vec);

void read_ws(std::ifstream& from);

bool read_token(std::ifstream& from, const std::string& block_del);
	
void init_tag_names_table();
	
bool read_tags(std::ifstream& from, 
	std::map<std::string, char>& read_names_table);	

void map_tags(std::ifstream& from, int verbose,
	const std::map<std::string, char>& read_names_table);
	
bool insert_name(const std::string& name, const char& value,
	std::map<std::string, char>& read_names_table);
	
void error(const std::string& msg);

void message(const std::string& msg, bool mode);

std::vector<std::string> tag_names_table(0);

const bool debug = false;

template<typename ET>
struct Input_type_selector {};

template<>
struct Input_type_selector<CGAL::Gmpq> {
	typedef CGAL::Gmpq Input_type;
};

template<>
struct Input_type_selector<GMP::Double> {
	typedef double Input_type;
};

template<>
struct Input_type_selector<CGAL::Gmpz> {
	typedef CGAL::Gmpz Input_type;
};

template<typename ET,
	typename Is_linear,
 	typename Is_symmetric,
	typename Has_no_inequalities>
bool doIt(int verbose, std::ifstream& from) {
	typedef typename Input_type_selector<ET>::Input_type Input_type;
	typedef Rep<ET, Input_type, Is_linear, Is_symmetric, Has_no_inequalities> Repr;
	Matrix<Input_type>  A;
	typename Repr::Row_type* row_types;
	Vector<Input_type>  b;
	Vector<Input_type>  c;
	Matrix<Input_type>  D;
	std::vector<int> rel;
	CGAL::QPE_solver< Repr >              solver;
	bool instance_read, sol_solver_feasible, sol_solver_optimal;	

	instance_read = read_instance<Input_type>(from, A, rel, b, c, D);
	if (instance_read) {
		init_row_types<Repr>(rel, row_types);
		solver.set_verbosity( verbose);
		solver.set( A.size(), rel.size(),
			typename Repr::Vector_iterator( A.begin(), Begin<Input_type>()),
			b.begin(), c.begin(),
			typename Repr::Vector_iterator( D.begin(), Begin<Input_type>()),
			row_types);
		CGAL::QPE_full_exact_pricing<Repr>  strategy;
		solver.set_pricing_strategy( strategy);
		solver.init();
		solver.solve();
		sol_solver_feasible = solver.is_solution_feasible();
		std::cout << "feasible: " << sol_solver_feasible << std::endl;
		sol_solver_optimal = solver.is_solution_optimal();
		std::cout << "optimal: " << sol_solver_optimal << std::endl;
		delete[] row_types;
		return true;
	} else {
		error("could not read problem instance");
		return false;
	}
}
 
 

int main( int argc, char** argv) {
	int verbose = 0;

	if ( argc > 2) {
		verbose = atoi( argv[ 1]);
	} else {
		return 1;
	}
	init_tag_names_table();
	std::map<std::string, char> read_names_table;
	
	for (int i = 2; i < argc; ++i) {
		std::ifstream from(argv[i]);
		if (!from) {
			std::cout << "could not open file: " << argv[i] << "\n";
		} else {
			std::cout << "processing file: " << argv[i] << "\n";
			read_tags(from, read_names_table);
			map_tags(from, verbose, read_names_table);
			from.close();
		}
	}

	
     /*
     Rep::Row_type*  row_types = new Rep::Row_type[ 1];
      row_types[ 0] = Rep::EQUAL;

      CGAL::QPE_solver<Rep>              solver;
      solver.set_verbosity( verbose);
      solver.set( 2, 1,
		  Vector_iterator( A.begin(), Begin()), b.begin(),
		  c.begin(), Vector_iterator( D.begin(), Begin()),
		  row_types);
      {
	CGAL::QPE_full_exact_pricing<Rep>  strategy;
	solver.set_pricing_strategy( strategy);
	solver.init();
	solver.solve();
	std::cerr << "-----------------------------------------------------\n";
      }
      {
	CGAL::QPE_partial_exact_pricing<Rep>  strategy;
	solver.set_pricing_strategy( strategy);
	solver.init();
	solver.solve();
	std::cerr << "-----------------------------------------------------\n";
      }
      {
	CGAL::QPE_full_filtered_pricing<Rep>  strategy;
	solver.set_pricing_strategy( strategy);
	solver.init();
	solver.solve();
	std::cerr << "-----------------------------------------------------\n";
      }
      {
	CGAL::QPE_partial_filtered_pricing<Rep>  strategy;
	solver.set_pricing_strategy( strategy);
	solver.init();
	solver.solve();
	std::cerr << "-----------------------------------------------------\n";
      }
      */
    
    return 0;
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
bool read_instance(std::ifstream& from, Matrix<T>& a,
	std::vector<int>& rel, Vector<T>& b,
	Vector<T>& c, Matrix<T>& d) {
	
	int m, n;
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
bool read_matrix(std::ifstream& from, int m, int n,
	//typename Vector<T>::size_type m,
	//typename Vector<T>::size_type n,
	Matrix<T>& mat) {
	
	std::ostringstream msgstr; 
	bool entry_read;	
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
	bool entry_read;
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
	bool entry_read;
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
	return entry_read;
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

void map_tags(std::ifstream& from, int verbose,
	const std::map<std::string, char>& read_names_table) {
	int offset;
	std::map<std::string, char>::const_iterator p;
	int index = 0;
	offset = tag_names_table.size() - 1;
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
				default :	error("impossible tag value");
				}
				
			} else { // 0,1
				switch (p->second) {
				case '0':	;
						break;
				case '1':	index += (1 << (offset - i));
						break;
				default :	error("impossible tag value");
				}
			}
		} else {
			error("unspecified tags");
			index = -1;
		}
	}
	switch (index) {
	/*
	case  0: 	doIt<CGAL::Gmpq,CGAL::Tag_false,CGAL::Tag_false,
				CGAL::Tag_false>(verbose, from);
			break;
	case  1:	doIt<CGAL::Gmpq,CGAL::Tag_false,CGAL::Tag_false,
				CGAL::Tag_true>(verbose, from);
			break;
	*/
	case  2:	doIt<CGAL::Gmpq,CGAL::Tag_false,CGAL::Tag_true,
				CGAL::Tag_false>(verbose, from);
			break;
	/*
	case  3:	doIt<CGAL::Gmpq,CGAL::Tag_false,CGAL::Tag_true,
				CGAL::Tag_true>(verbose, from);
			break;
	case  4:	doIt<CGAL::Gmpq,CGAL::Tag_true,CGAL::Tag_false,
				CGAL::Tag_false>(verbose, from);
			break;
	case  5:	doIt<CGAL::Gmpq,CGAL::Tag_true,CGAL::Tag_false,
				CGAL::Tag_true>(verbose, from);
			break;
	case  6:	doIt<CGAL::Gmpq,CGAL::Tag_true,CGAL::Tag_true,
				CGAL::Tag_false>(verbose, from);
			break;
	case  7:	doIt<CGAL::Gmpq,CGAL::Tag_true,CGAL::Tag_true,
				CGAL::Tag_true>(verbose, from);
			break;
	case  8:	doIt<GMP::Double,CGAL::Tag_false,CGAL::Tag_false,
				CGAL::Tag_false>(verbose, from);
			break;
	case  9:	doIt<GMP::Double,CGAL::Tag_false,CGAL::Tag_false,
				CGAL::Tag_true>(verbose, from);
			break;
	case 10:	doIt<GMP::Double,CGAL::Tag_false,CGAL::Tag_true,
				CGAL::Tag_false>(verbose, from);
			break;
	case 11:	doIt<GMP::Double,CGAL::Tag_false,CGAL::Tag_true,
				CGAL::Tag_true>(verbose, from);
			break;
	case 12:	doIt<GMP::Double,CGAL::Tag_true,CGAL::Tag_false,
				CGAL::Tag_false>(verbose, from);
			break;
	case 13:	doIt<GMP::Double,CGAL::Tag_true,CGAL::Tag_false,
				CGAL::Tag_true>(verbose, from);
			break;
	case 14:	doIt<GMP::Double,CGAL::Tag_true,CGAL::Tag_true,
				CGAL::Tag_false>(verbose, from);
			break;
	case 15:	doIt<GMP::Double,CGAL::Tag_true,CGAL::Tag_true,
				CGAL::Tag_true>(verbose, from);
			break;
	case 16:	doIt<CGAL::Gmpz,CGAL::Tag_false,CGAL::Tag_false,
				CGAL::Tag_false>(verbose, from);
			break;
	case 17:	doIt<CGAL::Gmpz,CGAL::Tag_false,CGAL::Tag_false,
				CGAL::Tag_true>(verbose, from);
			break;
	case 18:	doIt<CGAL::Gmpz,CGAL::Tag_false,CGAL::Tag_true,
				CGAL::Tag_false>(verbose, from);
			break;
	case 19:	doIt<CGAL::Gmpz,CGAL::Tag_false,CGAL::Tag_true,
				CGAL::Tag_true>(verbose, from);
			break;
	case 20:	doIt<CGAL::Gmpz,CGAL::Tag_true,CGAL::Tag_false,
				CGAL::Tag_false>(verbose, from);
			break;
	case 21:	doIt<CGAL::Gmpz,CGAL::Tag_true,CGAL::Tag_false,
				CGAL::Tag_true>(verbose, from);
			break;
	case 22:	doIt<CGAL::Gmpz,CGAL::Tag_true,CGAL::Tag_true,
				CGAL::Tag_false>(verbose, from);
			break;
	case 23:	doIt<CGAL::Gmpz,CGAL::Tag_true,CGAL::Tag_true,
				CGAL::Tag_true>(verbose, from);
			break;
	*/
	default:	;
	}
}

void init_tag_names_table() {
	tag_names_table.clear();
	tag_names_table.push_back("input_data_type");
	tag_names_table.push_back("is_linear");
	tag_names_table.push_back("is_symmetric");
	tag_names_table.push_back("has_no_inequalities");	
}


bool insert_name(const std::string& name, const char& value,
	std::map<std::string, char>& read_names_table) {
	
	std::pair<std::string, char> tmp(name, value);
	std::pair<std::map<std::string, char>::iterator,
		bool> p = read_names_table.insert(tmp);
	return p.second;
}

 
void read_ws(std::ifstream& from) {
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

// ===== EOF ==================================================================
