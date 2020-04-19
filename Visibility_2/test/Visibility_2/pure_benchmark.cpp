// Copyright (c) 2013 Technical University Braunschweig (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s):  Francisc Bungiu <fbungiu@gmail.com>
//             Michael Hemmer <michael.hemmer@cgal.org>
//             Kan Huang      <huangkandiy@gmail.com

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Rotational_sweep_visibility_2.h>
#include <CGAL/test_model_methods.h>
#include <CGAL/test_utils.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel     Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                    Traits_2;
typedef Traits_2::Point_2                                              Point_2;
typedef Traits_2::X_monotone_curve_2                                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                                  Arrangement_2;

template <class Visibility>
void deploy_pure_benchmark(const CGAL::Query_choice& qchoice, std::ifstream& input) {
  Visibility v;
  CGAL::pure_benchmark<Visibility>
      (v, qchoice, input);
}


template <class Regularization_category>
void benchmark_one_class(std::string name, const CGAL::Query_choice& qchoice, std::string input_arr_file) {
  std::ifstream input(input_arr_file.c_str());
  if (name == "S")
    deploy_pure_benchmark<CGAL::Simple_polygon_visibility_2<Arrangement_2, Regularization_category> > (qchoice, input);
  if (name == "T")
    deploy_pure_benchmark<CGAL::Triangular_expansion_visibility_2<Arrangement_2, Regularization_category> > (qchoice, input);
  if (name == "R")
    deploy_pure_benchmark<CGAL::Rotational_sweep_visibility_2<Arrangement_2, Regularization_category> > (qchoice, input);
}

void print_usage() {
  std::cout << "Usage: ./pure_benchmark [filename] [Class type] [Query type] [Regularize]\n";
  std::cout << "where [Class type] could be S(simple), R(rotational sweep), and T(triangular), indicating which class you want to test.\n";
  std::cout << "[Query type] can be: {vertex, edge, face}.\n";
  std::cout << "[Regularize] can be: {true, false}.\n";
}

int main(int argc, char* argv[]) {

  CGAL::Query_choice qchoice = CGAL::FACE;
  std::string regularization_tag("true");
  if (argc > 1) {
    std::string input_arr_file(argv[1]);

    if (argc == 2) {
      std::cout << "NAME TAG PreProTime NQueries TimeQueries TotalTime QAVE TAVE" << std::endl;
      benchmark_one_class<CGAL::Tag_true>(std::string("T"), CGAL::VERTEX, input_arr_file);
      benchmark_one_class<CGAL::Tag_true>(std::string("T"), CGAL::EDGE, input_arr_file);
      benchmark_one_class<CGAL::Tag_true>(std::string("T"), CGAL::FACE, input_arr_file);
      benchmark_one_class<CGAL::Tag_false>(std::string("T"), CGAL::VERTEX, input_arr_file);
      benchmark_one_class<CGAL::Tag_false>(std::string("T"), CGAL::EDGE, input_arr_file);
      benchmark_one_class<CGAL::Tag_false>(std::string("T"), CGAL::FACE, input_arr_file);
      std::cout << std::endl;
      benchmark_one_class<CGAL::Tag_true>(std::string("S"), CGAL::VERTEX, input_arr_file);
      benchmark_one_class<CGAL::Tag_true>(std::string("S"), CGAL::EDGE, input_arr_file);
      benchmark_one_class<CGAL::Tag_true>(std::string("S"), CGAL::FACE, input_arr_file);
      benchmark_one_class<CGAL::Tag_false>(std::string("S"), CGAL::VERTEX, input_arr_file);
      benchmark_one_class<CGAL::Tag_false>(std::string("S"), CGAL::EDGE, input_arr_file);
      benchmark_one_class<CGAL::Tag_false>(std::string("S"), CGAL::FACE, input_arr_file);
      std::cout << std::endl;
      benchmark_one_class<CGAL::Tag_true>(std::string("R"), CGAL::VERTEX, input_arr_file);
      benchmark_one_class<CGAL::Tag_true>(std::string("R"), CGAL::EDGE, input_arr_file);
      benchmark_one_class<CGAL::Tag_true>(std::string("R"), CGAL::FACE, input_arr_file);
      benchmark_one_class<CGAL::Tag_false>(std::string("R"), CGAL::VERTEX, input_arr_file);
      benchmark_one_class<CGAL::Tag_false>(std::string("R"), CGAL::EDGE, input_arr_file);
      benchmark_one_class<CGAL::Tag_false>(std::string("R"), CGAL::FACE, input_arr_file);
      std::cout << std::endl;
      benchmark_one_class<CGAL::Tag_true>(std::string("PR"), CGAL::VERTEX, input_arr_file);
      benchmark_one_class<CGAL::Tag_true>(std::string("PR"), CGAL::EDGE, input_arr_file);
      benchmark_one_class<CGAL::Tag_true>(std::string("PR"), CGAL::FACE, input_arr_file);
      benchmark_one_class<CGAL::Tag_false>(std::string("PR"), CGAL::VERTEX, input_arr_file);
      benchmark_one_class<CGAL::Tag_false>(std::string("PR"), CGAL::EDGE, input_arr_file);
      benchmark_one_class<CGAL::Tag_false>(std::string("PR"), CGAL::FACE, input_arr_file);
    } else if (argc == 5) {
      qchoice = CGAL::FACE;
      std::string query_type(argv[3]);
      if (query_type == "vertex")
        qchoice = CGAL::VERTEX;
      else if (query_type == "edge")
        qchoice = CGAL::EDGE;
      else if (query_type == "face")
        qchoice = CGAL::FACE;
      else {
        std::cout<<"query type is not matched.\n";
        return 0;
      }

      regularization_tag = argv[4];
      std::string classname(argv[2]);
      if (regularization_tag == "true") {
        benchmark_one_class<CGAL::Tag_true>(classname, qchoice, input_arr_file);
      }
      else {
        benchmark_one_class<CGAL::Tag_false>(classname, qchoice, input_arr_file);
      }

      return 0;
    }
    else {
      print_usage();
      exit(0);
    }
  }
  else {
    print_usage();
    exit(0);
  }
  return 0;
}
