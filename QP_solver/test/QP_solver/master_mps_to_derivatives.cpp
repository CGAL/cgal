// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
#include <boost/shared_ptr.hpp>

#include <CGAL/Exact_rational.h>
#include <CGAL/MP_Float.h>
#include <CGAL/QP_solver/QP_solver.h>
#include <CGAL/QP_solver/QP_full_exact_pricing.h>
#include <CGAL/QP_solver/QP_exact_bland_pricing.h>
#include <CGAL/QP_solver/QP_partial_exact_pricing.h>
#include <CGAL/QP_solver/QP_full_filtered_pricing.h>
#include <CGAL/QP_solver/QP_partial_filtered_pricing.h>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

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
#ifdef CGAL_USE_GMP  
  template<>
  struct MPS_type_name<CGAL::Gmpq> {
    static const char *name() { return "rational"; }
  };  
#endif
  
#ifdef CGAL_USE_LEDA
  template<>
  struct MPS_type_name<leda::rational> {
    static const char *name() { return "rational"; }
  };  
#endif
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

#ifdef CGAL_USE_GMP
  template<>
  struct IT_to_ET<int> {
    typedef CGAL::Gmpz ET;
  };
  
  template<>
  struct IT_to_ET<CGAL::Gmpq> {
    typedef CGAL::Gmpq ET;
  };
#endif

#ifdef CGAL_USE_LEDA
  template<>
  struct IT_to_ET<int> {
    typedef leda::integer ET;
  };
  
  template<>
  struct IT_to_ET<leda::rational> {
    typedef leda::rational ET;
  }; 
#endif
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
      <typename 
      std::iterator_traits<typename QP::U_iterator>::value_type>::name();
    if (tn != 0)
      out << "* Number-type: " << tn << "\n";
  } else
      out << "* Number-type: " << number_type << "\n";
  out << "* Description: " << description << "\n"
      << "* Generated-by: " << generator_name << "\n";

  // print in qmatrix format
  CGAL::print_quadratic_program(out, qp, problem_name);
}

boost::shared_ptr<std::ofstream>
create_output_file(const char *filename, // Note: "Bernd3" and not
					 // "Bernd3.mps".
		   const char *directory,
		   const char *suffix)   // Note: "shifted"
					 // and not "_shifted".
{
  std::string new_name = std::string(directory) +
    std::string("/") + std::string(filename) + std::string("_") +
    std::string(suffix) + std::string(".mps");
  return boost::shared_ptr<std::ofstream>(new std::ofstream(new_name.c_str(),
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
void create_shifted_instance(const CGAL::Quadratic_program_from_mps <IT>& qp,
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
  const int n = qp.get_n();
  const int m = qp.get_m();

  // offset vector:
  std::vector<IT> v(n);
  for (int i=0; i<n; ++i)
    v[i] = i+1;

  // compute A v into Av;
  std::vector<IT> Av(m, IT(0));
  for (int i=0; i<m; ++i) 
    for (int j=0; j<n; ++j)
      Av[i] += (*(qp.get_a()+j))[i] * v[j];

  // compute - 2 v^T D into mvTD:
  std::vector<IT> mvTD(n, IT(0));  // -2D^Tv
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j)
      mvTD[i] 
	+= ( j <= i ? (*(qp.get_d()+i))[j] : (*(qp.get_d()+j))[i]) * v[j];
    mvTD[i] *= -1;
  }

  // output:
  using boost::make_transform_iterator;
  using boost::make_zip_iterator;
  boost::shared_ptr<std::ofstream> out = create_output_file(file, dir, "shifted");

  write_MPS(*out,
		  "", // deduce number-type
		  "Shifted instance of original file",
		  "master_mps_to_derivatives-create_shifted_instance",
		  qp.get_problem_name(), 
		  CGAL::make_quadratic_program_from_iterators(
     n, 
     m, 
     qp.get_a(), 
     make_transform_iterator(
			     make_zip_iterator(boost::make_tuple(qp.get_b(),Av.begin())),
			     tuple_add<IT>()),
     qp.get_r(), 
     qp.get_fl(),
     make_transform_iterator(
			      make_zip_iterator(boost::make_tuple(qp.get_l(),v.begin())),
			      tuple_add<IT>()),
     qp.get_fu(),
     make_transform_iterator(
			     make_zip_iterator(boost::make_tuple(qp.get_u(),v.begin())),
			     tuple_add<IT>()),
     qp.get_d(),
     make_transform_iterator(
			   make_zip_iterator(boost::make_tuple(qp.get_c(),mvTD.begin())),
			   tuple_add<IT>()), 
     qp.get_c0()
     )
);
  out->close();
}

template<typename IT,   // input number type
	 typename ET>   // exact number type compatible with IT (ET is used,
                        // for instance, by QP in the query methods
                        // has_equalities_only_and_full_rank())
void create_free_instance(CGAL::Quadratic_program_from_mps<IT>& qp,
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
  const unsigned int n = qp.get_n();
  unsigned int m = qp.get_m();

  // add rows to A and corresponding entries to b:
  for (unsigned int i=0; i<n; ++i) {
    if (*(qp.get_fl()+i)) {                        // x >= l
      // add a row to A:
      for (unsigned int j=0; j<n; ++j)
	qp.set_a (j, m, (i==j? 1 : 0)); 

      // add corresponding entry to b:
      qp.set_b(m, qp.get_l()[i]);

      // add corresponding row type:
      qp.set_r(m, CGAL::LARGER);
      ++m;
    }
    qp.set_l(i, false);                           // variable becomes free
    if (*(qp.get_fu()+i)) {                        // x <= u
      // add a row to A:
      for (unsigned int j=0; j<n; ++j)
	qp.set_a (j, m ,(i==j? 1 : 0)); 

      // add corresponding entry to b:
      qp.set_b(m, qp.get_u()[i]);

      // add corresponding row type:
      qp.set_r(m, CGAL::SMALLER);
      ++m;
    }
    qp.set_u(i, false);                         // variable becomes free
  }
  // output:
  boost::shared_ptr<std::ofstream> out = create_output_file(file, dir, "free");
  write_MPS(*out,
		  "", // deduce number-type
		  "Freed instance of original file",
		  "master_mps_to_derivatives-create_free_instance",
		  qp.get_problem_name(),
	    qp);
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
  typedef typename QP_from_mps_detail::IT_to_ET<IT>::ET ET;
  typedef CGAL::Quadratic_program_from_mps<IT> QP;
  QP qp(f);

  // check for format errors in MPS file:
  if (!qp.is_valid()) {
    msg = "Input is not a valid MPS file: " + qp.get_error();
    return false;
  }
  cerr << "    MPS-file successfully input.\n";

  // no derivatives if comment says so
  if (qp.get_comment().find(std::string("Derivatives: none"))!=std::string::npos) 
    cerr << "    No derivatives made.\n";
  else {
    // derivates:
    create_shifted_instance<IT, ET>(qp, file, dir);
    create_free_instance<IT, ET>(qp, file, dir);
    // Note: insert additional derivative routines here! Your routine may use
    // create_output_file() to create the output file.
  }
  // cleanup:
  f.close();
  return true;
}

int main(const int argnr, const char **argv) {
  typedef CGAL::Exact_rational Rational;

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
