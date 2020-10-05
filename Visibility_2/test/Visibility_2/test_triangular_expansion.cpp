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


#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/test_model_methods.h>
#include <CGAL/test_utils.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Arr_extended_dcel.h>

#include <iostream>
#include <fstream>

int main() {
  {
    typedef CGAL::Exact_rational                                  Number_type;
    typedef CGAL::Cartesian<Number_type>                          Kernel;
    typedef CGAL::Arr_segment_traits_2<Kernel>                    Traits_2;
    typedef CGAL::Arrangement_2<Traits_2>                         Arrangement_2;
    {
      typedef CGAL::Triangular_expansion_visibility_2<
                    Arrangement_2 , CGAL::Tag_true>               Visibility_2;
      CGAL::test_model_methods<Visibility_2,Arrangement_2>();
      CGAL::run_tests<Visibility_2,Arrangement_2>(22, 2);
    }
    {
      typedef CGAL::Triangular_expansion_visibility_2<
                    Arrangement_2 ,  CGAL::Tag_false>             Visibility_2;
      CGAL::test_model_methods<Visibility_2,Arrangement_2>();
      CGAL::run_tests<Visibility_2,Arrangement_2>(22, 2);
    }
  }
  {
    typedef CGAL::Exact_predicates_exact_constructions_kernel     Kernel;
    typedef CGAL::Arr_segment_traits_2<Kernel>                    Traits_2;
    typedef CGAL::Arrangement_2<Traits_2>                         Arrangement_2;
    {
      typedef CGAL::Triangular_expansion_visibility_2<
                    Arrangement_2 , CGAL::Tag_true>               Visibility_2;
      CGAL::test_model_methods<Visibility_2,Arrangement_2>();
      CGAL::run_tests<Visibility_2,Arrangement_2>(22, 2);
    }{
      typedef CGAL::Triangular_expansion_visibility_2<
                    Arrangement_2 , CGAL::Tag_false>              Visibility_2;
      CGAL::test_model_methods<Visibility_2,Arrangement_2>();
      CGAL::run_tests<Visibility_2,Arrangement_2>(22, 2);
    }
  }

  {
    // test Visibility_arrangement_type with extended DCEL
    typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
    typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
    typedef CGAL::Arrangement_2<Traits_2>                           ARR;
    typedef CGAL::Arr_extended_dcel<Traits_2, bool, bool, bool>     EDCEL;
    typedef CGAL::Arrangement_2<Traits_2, EDCEL>                    EARR;
    {
      typedef CGAL::Triangular_expansion_visibility_2<ARR,CGAL::Tag_true> Visibility_2;
      CGAL::test_model_methods<Visibility_2,EARR>();
      CGAL::run_tests<Visibility_2,EARR>(22, 2);
    }{
      typedef CGAL::Triangular_expansion_visibility_2<ARR,CGAL::Tag_false> Visibility_2;
      CGAL::test_model_methods<Visibility_2,EARR>();
      CGAL::run_tests<Visibility_2,EARR>(22, 2);
    }
  }

  return 0;
}
