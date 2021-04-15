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
//             Kan Huang      <huangkandiy@gmail.com>
//             Ning Xu        <longyin0904@gmail.com>


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

template <class Visibility_fst, class Visibility_snd>
void deploy_benchmark(CGAL::Query_choice& qchoice, std::ifstream& input) {
  Visibility_fst v1;
  Visibility_snd v2;
  CGAL::simple_benchmark<Visibility_fst, Visibility_snd>
      (v1, v2, qchoice, input);
}

template <class Visibility_fst, class Regularization_category>
void define_snd_class(std::string name2, CGAL::Query_choice& qchoice, std::ifstream& input){
  if (name2 == "S")
    deploy_benchmark<Visibility_fst, CGAL::Simple_polygon_visibility_2<Arrangement_2, Regularization_category> >
        (qchoice, input);
  if (name2 == "T")
    deploy_benchmark<Visibility_fst, CGAL::Triangular_expansion_visibility_2<Arrangement_2, Regularization_category> >
        (qchoice, input);
  if (name2 == "R")
    deploy_benchmark<Visibility_fst, CGAL::Rotational_sweep_visibility_2<Arrangement_2, Regularization_category> >
        (qchoice, input);
}

template <class Regularization_category>
void benchmark_two_classes(std::string name1, std::string name2, CGAL::Query_choice& qchoice, std::ifstream& input) {
  if (name1 == "S")
    define_snd_class<CGAL::Simple_polygon_visibility_2<Arrangement_2, Regularization_category>, Regularization_category> (name2, qchoice, input);
  if (name1 == "T")
    define_snd_class<CGAL::Triangular_expansion_visibility_2<Arrangement_2, Regularization_category>, Regularization_category> (name2, qchoice, input);
  if (name1 == "R")
    define_snd_class<CGAL::Rotational_sweep_visibility_2<Arrangement_2, Regularization_category>, Regularization_category> (name2, qchoice, input);
}

void print_usage() {
  std::cout << "Usage: ./simple_benchmark [filename] [Class type 1] [Class type 2] [Query type] [Regularize]\n";
  std::cout << "where [Class type] could be S(simple), R(rotational sweep) and T(triangular), indicating which classes you want to test.\n";
  std::cout << "[Query type] can be: {vertex, edge, face}.\n";
  std::cout << "[Regularize] can be: {true, false}.\n";
}

int main(int argc, char* argv[]) {

  CGAL::Query_choice qchoice = CGAL::FACE;
  std::string regularization_tag("true");
  if (argc > 1) {
    std::string input_arr_file(argv[1]);
    std::ifstream input(input_arr_file.c_str());
    if (argc == 6) {
      qchoice = CGAL::FACE;
      std::string query_type(argv[4]);
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

      regularization_tag = argv[5];
      std::string class1(argv[2]), class2(argv[3]);
      if (regularization_tag == "true") {
        benchmark_two_classes<CGAL::Tag_true>(class1, class2, qchoice, input);
      }
      else {
        benchmark_two_classes<CGAL::Tag_false>(class1, class2, qchoice, input);
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
