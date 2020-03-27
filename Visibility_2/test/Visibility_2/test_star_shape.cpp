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
// Author(s):  Kan Huang <huangkandiy@gmail.com>
//             Francisc Bungiu <fbungiu@gmail.com>
//             Michael Hemmer <michael.hemmer@cgal.org>


#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/test_model_methods.h>
#include <CGAL/test_utils.h>


#include <iostream>
#include <fstream>

int main(int argc, char* argv[]) {
{
  typedef CGAL::Exact_rational                                  Number_type;
  typedef CGAL::Cartesian<Number_type>                                                         Kernel;
  typedef CGAL::Arr_segment_traits_2<Kernel>                    Traits_2;
  typedef CGAL::Arrangement_2<Traits_2>                                                        Arrangement_2;
  typedef CGAL::Simple_polygon_visibility_2<
                Arrangement_2, CGAL::Tag_false>     Simple_polygon_visibility_2;

  typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2,  CGAL::Tag_false>
    Triangular_expansion_visibility_2;

  if (argc == 2) {
    Simple_polygon_visibility_2 simple_visibility;
    const CGAL::Query_choice qchoice = CGAL::FACE;
    std::string input_arr_file(argv[1]);
    std::ifstream input(input_arr_file.c_str());
    CGAL::test_star_shape<Simple_polygon_visibility_2>
                  (simple_visibility, qchoice, input);
    return 0;
  }

  if (argc == 3) {
    const CGAL::Query_choice qchoice = CGAL::FACE;
    std::string input_arr_file(argv[1]);
    std::ifstream input(input_arr_file.c_str());
    std::string class_name(argv[2]);
    if ( class_name == "simple") {
      Simple_polygon_visibility_2 simple_visibility;
      CGAL::test_star_shape<Simple_polygon_visibility_2>
                    (simple_visibility, qchoice, input);
      return 0;
    }
    if (class_name == "triangular") {
      Triangular_expansion_visibility_2 triangular_visibility;
      CGAL::test_star_shape<Triangular_expansion_visibility_2>
                    (triangular_visibility, qchoice, input);
      return 0;
    }
    std::cout<<"no type is matched.\n";
    return 0;
  }

  if (argc == 4) {
    std::string input_arr_file(argv[1]);
    std::ifstream input(input_arr_file.c_str());
    CGAL::Query_choice qchoice;
    std::string query_type(argv[3]);
    if (query_type == "vertex")
      qchoice = CGAL::VERTEX;
    else {
      if (query_type == "edge")
        qchoice = CGAL::EDGE;
      else {
        if (query_type == "face")
          qchoice = CGAL::FACE;
        else {
          std::cout<<"query type is not matched.\n";
          return 0;
        }
      }
    }
    std::string class_name(argv[2]);
    if (class_name == "simple") {
      Simple_polygon_visibility_2 simple_visibility;
      CGAL::test_star_shape<Simple_polygon_visibility_2>
                    (simple_visibility, qchoice, input);
      return 0;
    }
    if (class_name == "triangular") {
      Triangular_expansion_visibility_2 triangular_visibility;
      CGAL::test_star_shape<Triangular_expansion_visibility_2>
                    (triangular_visibility, qchoice, input);
      return 0;
    }
    std::cout<<"no type is matched.\n";
    return 0;

  }


  std::cout << "Usage: ./test_star_shape [filename] [Class type] [Query type]\n";
  std::cout << "where [Class type] could be simple, naive and triangular, indicating which class you want to test.\n";
  std::cout << "[Query type] could be vertex, edge, face.\n";
  std::cout << "The default value of [Query type] is face. The default value of [Class type] is simple.\n";
  exit(0);
}

}
