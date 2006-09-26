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
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <CGAL/Gmpq.h>
#include <CGAL/QP_solver.h>
#include <CGAL/QP_solver/QP_full_exact_pricing.h>
#include <CGAL/QP_solver/QP_exact_bland_pricing.h>
#include <CGAL/QP_solver/QP_partial_exact_pricing.h>
#include <CGAL/QP_solver/QP_full_filtered_pricing.h>
#include <CGAL/QP_solver/QP_partial_filtered_pricing.h>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

using CGAL::Tag_true;
using CGAL::Tag_false;

typedef Tag_false Sparse_D;
typedef Tag_false Sparse_A;
typedef Tag_false Is_linear;

// Routines to output to MPS format:
namespace QP_from_mps_detail {

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
  struct MPS_type_name<CGAL::Gmpq> {
    static const char *name() { return "rational"; }
  };  

  template<>
  struct MPS_type_name<CGAL::Quotient<CGAL::MP_Float> > {
    static const char *name() { return "rational"; }
  }; template<typename IT>
  struct IT_to_ET {
  };
  
  template<>
  struct IT_to_ET<double> {
    typedef CGAL::MP_Float ET;
  };
  
  template<>
  struct IT_to_ET<int> {
    typedef CGAL::Gmpz ET;
  };
  
  template<>
  struct IT_to_ET<CGAL::Gmpq> {
    typedef CGAL::Gmpq ET;
  }; 

  template<>
  struct IT_to_ET<CGAL::Quotient<CGAL::MP_Float> > {
    typedef CGAL::Quotient<CGAL::MP_Float> ET;
  };

} // QP_from_mps_detail

template<typename QP>
void write_MPS(std::ostream& out,
	       const std::string& number_type, // pass "" to deduce
					       // the number-type from 
                                               // U_iterator::value_type
	       const std::string& description,
	       const std::string& generator_name,
	       const std::string& problem_name,
	       const QP& qp)
{
  // output header:
  if (number_type.length() == 0) {
    const char *tn = QP_from_mps_detail::MPS_type_name
      <typename QP::U_iterator::value_type>::name();
    if (tn != 0)
      out << "* Number-type: " << tn << "\n";
  } else
      out << "* Number-type: " << number_type << "\n";
  out << "* Description: " << description << "\n"
      << "* Generated-by: " << generator_name << "\n";

  CGAL::print_quadratic_program(out, qp, problem_name);
}

std::auto_ptr<std::ofstream>
create_output_file(const char *filename, // Note: "Bernd3" and not
					 // "Bernd3.mps".
		   const char *directory,
		   const char *suffix)   // Note: "shifted"
					 // and not "_shifted".
{
  std::string new_name = std::string(directory) +
    std::string("/") + std::string(filename) + std::string("_") +
    std::string(suffix) + std::string(".mps");
  return std::auto_ptr<std::ofstream>(new std::ofstream(new_name.c_str(),
							std::ios_base::trunc |
							std::ios_base::out));
}

template<typename NT>
struct tuple_add : 
  public std::unary_function<const boost::tuple<NT, NT>&, NT>
{
  NT operator()(const boost::tuple<NT, NT>& t) const
  {
    return boost::tuples::get<0>(t) + boost::tuples::get<1>(t);
  }
};

template<typename IT,   // input number type
	 typename ET>   // exact number type compatible with IT (ET is used,
                        // for instance, by QP in the query methods
                        // has_equalities_only_and_full_rank())
void create_shifted_instance(CGAL::QP_from_mps
			     <IT, Is_linear, Sparse_D, Sparse_A>& qp,
			     const char *path,
			     const char *file,   // Note: "Bernd3" and
					         // not "Bernd3.mps".
			     const char *dir)
{
  // This routine implements the following transformation:
  //
  //   A  -> A
  //   b  -> b + A v
  //   c  -> c - 2 v^T 
  //   D  -> D
  //   row_types -> row_types
  //   l  -> l + v
  //   u  -> u + v
  //   fl -> fl
  //   fu -> fu
  //
  // where v = [1,...,n]^T.

  // extract data from qp:
  const int n = qp.n();
  const int m = qp.m();

  // offset vector:
  std::vector<IT> v(n);
  for (int i=0; i<n; ++i)
    v[i] = i+1;

  // compute A v into Av;
  std::vector<IT> Av(m, IT(0));
  for (int i=0; i<m; ++i) 
    for (int j=0; j<n; ++j)
      Av[i] += qp.a()[j][i] * v[j];

  // compute - 2 v^T D into mvTD:
  std::vector<IT> mvTD(n, IT(0));  // -2D^Tv
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j)
      mvTD[i] += qp.d()[j][i] * v[j];
    mvTD[i] *= -2;
  }

  // output:
  using boost::make_transform_iterator;
  using boost::make_zip_iterator;
  using boost::make_tuple;
  std::auto_ptr<std::ofstream> out = create_output_file(file, dir, "shifted");

  write_MPS(*out,
		  "", // deduce number-type
		  "Shifted instance of original file",
		  "master_mps_to_derivatives-create_shifted_instance",
		  qp.problem_name(), 
		  CGAL::make_QP_from_iterators(
     n, 
     m, 
     qp.a(), 
     make_transform_iterator(
			     make_zip_iterator(make_tuple(qp.b(),Av.begin())),
			     tuple_add<IT>()),
     qp.r(), 
     qp.fl(),
     make_transform_iterator(
			      make_zip_iterator(make_tuple(qp.l(),v.begin())),
			      tuple_add<IT>()),
     qp.fu(),
     make_transform_iterator(
			     make_zip_iterator(make_tuple(qp.u(),v.begin())),
			     tuple_add<IT>()),
     qp.d(),
     make_transform_iterator(
			   make_zip_iterator(make_tuple(qp.c(),mvTD.begin())),
			   tuple_add<IT>()), 
     qp.c0()
     )
);
// 		  n, m,
// 		  qp.a(),
// 		  make_transform_iterator(
// 		    make_zip_iterator(make_tuple(qp.b(),Av.begin())),
// 		    tuple_add<IT>()),
// 		  make_transform_iterator(
// 		    make_zip_iterator(make_tuple(qp.c(),mvTD.begin())),
// 		    tuple_add<IT>()), qp.c0(),
// 		  qp.d(), qp.fu(), qp.fl(),
// 		  make_transform_iterator(
// 		    make_zip_iterator(make_tuple(qp.u(),v.begin())),
// 		    tuple_add<IT>()),
// 		  make_transform_iterator(
// 		    make_zip_iterator(make_tuple(qp.l(),v.begin())),
// 		    tuple_add<IT>()),
// 		  qp.r());
  out->close();
}

template<typename IT,   // input number type
	 typename ET>   // exact number type compatible with IT (ET is used,
                        // for instance, by QP in the query methods
                        // has_equalities_only_and_full_rank())
void create_free_instance(CGAL::QP_from_mps<IT, Is_linear,
			  Sparse_D, Sparse_A>& qp_,
			  const char *path,
			  const char *file,   // Note: "Bernd3" and
			                      // not "Bernd3.mps".
			  const char *dir)
{
  // This routine converts the given instance into an equivalent
  // problem where all bounds are modelled by additional rows of A and
  // where all variables are free.
  //
  // That is, the quantities c and D do not change, but A, b, and
  // row_types are augmented by at most 2n additional rows/entries
  // (and fl and fu are adjusted as well).

  // extract data from qp:
  const unsigned int n = qp_.n();
  const unsigned int m = qp_.m();

  // allocate storage (admittedly, I don't care about efficiency and
  // elegance here...):
  typedef CGAL::QP_from_mps<IT, Is_linear, Sparse_D, Sparse_A>    QP_MPS;
  typedef typename QP_MPS::Vector            Vector;
  typedef typename QP_MPS::B_iterator   Vector_iterator;
  typedef typename QP_MPS::A_Matrix          Matrix;
  typedef typename QP_MPS::A_Beginner        A_Beginner;
  typedef typename QP_MPS::A_iterator        A_iterator;
  typedef typename CGAL::Comparison_result Row_type;
  typedef typename QP_MPS::R_vector R_vector;

  // copy the qp
  QP_MPS qp (qp_);

  // copy some of its vectors (they get manipulated)
  Vector b           = qp.b_vector();
  R_vector row_types = qp.r_vector();

  // add rows to A and corresponding entries to b:
  int nr_of_rows_added = 0;
  for (unsigned int i=0; i<n; ++i) {
    if (*(qp.fl()+i)) {                        // x >= l
      // add a row to A:
      for (unsigned int j=0; j<n; ++j)
	qp.add_entry_in_A (j, b.size(),(i==j? 1 : 0)); 
        //A[j].push_back(i==j? 1 : 0);

      // add corresponding entry to b:
      b.push_back(qp.l()[i]);

      // add corresponding row type:
      row_types.push_back(CGAL::LARGER);

      ++nr_of_rows_added;
    }
    if (*(qp.fu()+i)) {                        // x <= u
      // add a row to A:
      for (unsigned int j=0; j<n; ++j)
	qp.add_entry_in_A (j, b.size(),(i==j? 1 : 0)); 
        //A[j].push_back(i==j? 1 : 0);

      // add corresponding entry to b:
      b.push_back(qp.u()[i]);

      // add corresponding row type:
      row_types.push_back(CGAL::SMALLER);

      ++nr_of_rows_added;
    }
  }

  // output:
  std::auto_ptr<std::ofstream> out = create_output_file(file, dir, "free");
  write_MPS(*out,
		  "", // deduce number-type
		  "Freed instance of original file",
		  "master_mps_to_derivatives-create_free_instance",
		  qp.problem_name(),
		  CGAL::make_QP_from_iterators (
		  n, m+nr_of_rows_added,
		  A_iterator(qp.A_matrix().begin(),A_Beginner()),
		  b.begin(),
		  row_types.begin(),
		  CGAL::Const_oneset_iterator<bool>(false), // fl
		  qp.l(),  // dummy
		  CGAL::Const_oneset_iterator<bool>(false), // fu
		  qp.u(),  // dummy
		  qp.d(),
		  qp.c(), qp.c0()));
  out->close();
}

template<typename IT>
bool create_derivatives(const char *path,
			const char *file,
			const char *dir, 
			std::string& msg)
{
  using std::cerr;
  using CGAL::Tag_true;
  using CGAL::Tag_false;

  // diagnostics:
  cerr << "  Trying to load input MPS using "
       << QP_from_mps_detail::MPS_type_name<IT>::name()
       << " number-type...\n";

  // open input file:
  std::ifstream f(path);
  if (!f) {
    cerr << "    Could not open file '" << path << "'.\n";
    return false;
  }

  // load QP instance:
  const int verbosity = 5;
  typedef typename QP_from_mps_detail::IT_to_ET<IT>::ET ET;
  typedef CGAL::QP_from_mps<IT, Is_linear, Sparse_D, Sparse_A> QP;
  QP qp(f,true,verbosity);

  // check for format errors in MPS file:
  if (!qp.is_valid()) {
    msg = "Input is not a valid MPS file: " + qp.error();
    return false;
  }
  cerr << "    MPS-file successfully input.\n";

  // no derivatives if comment says so
  if (qp.comment().find(std::string("Derivatives: none"))!=std::string::npos) 
    cerr << "    No derivatives made.\n";
  else {
    // derivates:
    create_shifted_instance<IT, ET>(qp, path, file, dir);
    create_free_instance<IT, ET>(qp, path, file, dir);
    // Note: insert additional derivative routines here! Your routine may use
    // create_output_file() to create the output file.
  }
  // cleanup:
  f.close();
  return true;
}

int main(const int argnr, const char **argv) {
  typedef CGAL::Gmpq Rational;

  // output usage information:
  if (argnr != 4) {
    std::cerr << "Usage: " << argv[0] << " path-to-master.mps name "
	      << "dest-directory\n\n"
	      << "Given a master MPS-file, this program constructs from it "
	      << "several other, derived\nMPS-files and saves them in the "
	      << "destination directory.\n\n"
	      << "The argument 'name' should be the basename (without "
	      << "extension) of the\nargument 'path-to-master.mps'.\n\n"
	      << "Usually, you do not have to call this program directly; "
	      << "./create_testsuite does\nit for you automagically.\n";
    return 0; // Note: 0 because otherwise the testsuite will fail...
  }

  // extract arguments:
  const char *path = argv[1];
  const char *file = argv[2];
  const char *dir  = argv[3];

  // As we do not know the number-type used in the MPS-file, we simply try
  // to load the MPS-file once with a double type, once with a Rational
  // type, and once with an int type.
  // we try the most special one first, so the order is 
  // int -> double -> rational
  std::string message;
  if (!create_derivatives<int>(path, file, dir, message))
    if (!create_derivatives<double>(path, file, dir, message))
      if (!create_derivatives<Rational>(path, file, dir, message)) {
	// Here, the MPS-file must be ill-formatted.
	std::cerr << "  " << message << "\n";
	return 2;
      }

  return 0;
}
