// Copyright (c) 2017 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_CREATION_OF_TEST_CASES_FOR_PATHS_H
#define CGAL_CREATION_OF_TEST_CASES_FOR_PATHS_H 1

#include<CGAL/Path_generators.h>

namespace CGAL {

template<typename Path>
void generate_one_positive_spur(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5]));
  extend_straight_positive(p, 6);
  extend_uturn_positive(p);
  extend_uturn_half_turn(p);
  extend_straight_positive(p, 4);
}

template<typename Path>
void generate_one_negative_spur(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5]));
  extend_straight_negative(p, 6);
  extend_uturn_negative(p);
  extend_uturn_half_turn(p);
  extend_straight_negative(p, 4);
}

template<typename Path>
void generate_cyclic_spur(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5]));
  extend_uturn_half_turn(p);
}

template<typename Path>
void generate_one_positive_bracket(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5]));
  extend_straight_positive(p, 3);
  extend_uturn_positive(p, 3);
  extend_uturn_positive(p);
  extend_straight_positive(p, 6);
  extend_uturn_positive(p);
  extend_uturn_positive(p, 3);
  extend_straight_positive(p, 2);
}

template<typename Path>
void generate_one_negative_bracket(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5]));
  extend_straight_negative(p, 3);
  extend_uturn_negative(p);
  extend_straight_negative(p, 6);
  extend_uturn_negative(p);
  extend_straight_negative(p, 2);
}

template<typename Path>
void generate_bracket_special1(Path& p, bool reverse)
{ // Case (x, 1, 2^r, 1)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[79]));
  extend_uturn_positive(p, 1);
  extend_straight_positive(p, 7);
  extend_uturn_positive(p);

  if (reverse)
  { p.reverse(); }
}

template<typename Path>
void generate_positive_bracket_special1(Path& p)
{ generate_bracket_special1(p, false); }

template<typename Path>
void generate_negative_bracket_special1(Path& p)
{ generate_bracket_special1(p, true); }

template<typename Path>
void generate_bracket_special2(Path& p, bool reverse)
{ // Case (1, 2^r)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[0]));
  extend_uturn_positive(p, 1);
  extend_straight_positive(p, 10);

  if (reverse)
  { p.reverse(); }
}

template<typename Path>
void generate_positive_bracket_special2(Path& p)
{ generate_bracket_special2(p, false); }

template<typename Path>
void generate_negative_bracket_special2(Path& p)
{ generate_bracket_special2(p, true); }

template<typename Path>
void generate_one_l_shape(Path& p)
{ // Generic case (... x -2^s -1 -2^t y ... ): here  (-2^2 -3 -2^8 -1 -2^5 -3 -2^3)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[0]));
  extend_straight_negative(p, 2);
  extend_uturn_negative(p, 3);
  extend_straight_negative(p, 8);
  extend_uturn_negative(p);
  extend_straight_negative(p, 5);
  extend_uturn_negative(p, 3);
  extend_straight_negative(p, 3);
}

template<typename Path>
void generate_l_shape_case2(Path& p)
{ // (... x -1 -2^t y ...): here (-2^2 -3 -1 -2^5 -3 -2^3)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[0]));
  extend_straight_negative(p, 2);
  extend_uturn_negative(p, 3);
  extend_uturn_negative(p);
  extend_straight_negative(p, 5);
  extend_uturn_negative(p, 3);
  extend_straight_negative(p, 3);
}

template<typename Path>
void generate_l_shape_case3(Path& p)
{ // (... x -2^s -1 y ...): here (-2^2 -3 -2^5 -1 -3 -2^3)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[0]));
  extend_straight_negative(p, 2);
  extend_uturn_negative(p, 3);
  extend_straight_negative(p, 5);
  extend_uturn_negative(p, 1);
  extend_uturn_negative(p, 3);
  extend_straight_negative(p, 3);
}

template<typename Path>
void generate_l_shape_case4(Path& p)
{ // (x -2^s -1 -2^t): here (-2^7 -1 -2^3 -4)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[35]));
  extend_uturn_negative(p, 4);
  extend_straight_negative(p, 7);
  extend_uturn_negative(p, 1);
  extend_straight_negative(p, 2);
}

template<typename Path>
void generate_l_shape_case5(Path& p)
{ // (x -1 -2^t): here (-4 -1 -2^12)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[27]));
  extend_uturn_negative(p, 1);
  extend_straight_negative(p, 12);
}

template<typename Path>
void generate_l_shape_case6(Path& p)
{ // (x -2^t -1): here (-4 -2^12 -1)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[13]));
  extend_straight_negative(p, 12);
  extend_uturn_negative(p, 1);
}

template<typename Path>
void generate_l_shape_case7(Path& p)
{ // (-3 -2^s -1 -2^t): here (-2^7 -1 -2^3 -3)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[194]));
  extend_straight_negative(p, 7);
  extend_uturn_negative(p, 1);
  extend_straight_negative(p, 3);
}

template<typename Path>
void generate_l_shape_case8(Path& p)
{ // (-2^l): here (-2^20)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[4]));
  extend_straight_negative(p, 19);
}

template<typename Path>
void generate_g1_v0_torus(Path& p)
{ // 1st generator
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[0]));
  extend_straight_negative(p, 4);
}

template<typename Path>
void generate_g1_v1_torus(Path& p)
{ // 1st generator v2
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[20]));
  extend_straight_negative(p, 4);
}

template<typename Path>
void generate_g1_v2_torus(Path& p)
{ // 1st generator v3
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[20]));
  extend_uturn_negative(p, 1);
  extend_straight_negative(p, 1);
  extend_uturn_positive(p, 1);
  extend_straight_negative(p, 2);
  extend_uturn_positive(p, 1);
  extend_straight_negative(p, 2);
  extend_uturn_negative(p, 1);
  extend_uturn_negative(p, 1);
}

template<typename Path>
void generate_g1_torus(Path& p, std::size_t i)
{
  assert(i<3);
  switch(i)
  {
    case 0: generate_g1_v0_torus(p); break;
    case 1: generate_g1_v1_torus(p); break;
    case 2: generate_g1_v2_torus(p); break;
    default: assert(false);
  }
}

template<typename Path>
void generate_null_cycle_v0_torus(Path& p)
{ // Empty cycle
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[20]));
  extend_uturn_positive(p, 1);
  extend_uturn_positive(p, 1);
  extend_uturn_positive(p, 1);
}

template<typename Path>
void generate_null_cycle_v1_torus(Path& p)
{ // Empty cycle v2
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[20]));
  extend_uturn_positive(p, 1);
  extend_straight_positive(p, 1);
  extend_uturn_positive(p, 1);
  extend_uturn_positive(p, 1);
  extend_straight_positive(p, 1);
}

template<typename Path>
void generate_null_cycle_v2_torus(Path& p)
{ // Empty cycle v3
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[20]));
  extend_uturn_negative(p, 1);
  extend_straight_negative(p, 3);
  extend_uturn_negative(p, 1);
  extend_straight_negative(p, 2);
  extend_uturn_negative(p, 1);
  extend_straight_negative(p, 3);
  extend_uturn_negative(p, 1);
  extend_straight_negative(p, 1);
}

template<typename Path>
void generate_null_cycle_torus(Path& p, std::size_t i)
{
  assert(i<3);
  switch(i)
  {
    case 0: generate_null_cycle_v0_torus(p); break;
    case 1: generate_null_cycle_v1_torus(p); break;
    case 2: generate_null_cycle_v2_torus(p); break;
    default: assert(false);
  }
}

template<typename Path>
void generate_g2_v0_torus(Path& p)
{ // 2nd generator v1
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[1]));
  extend_straight_positive(p, 4);
}

template<typename Path>
void generate_g2_v1_torus(Path& p)
{ // 2nd generator v2
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[13]));
  extend_straight_negative(p, 4);
}

template<typename Path>
void generate_g2_v2_torus(Path& p)
{ // 2nd generator v3
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[0]));
  extend_uturn_negative(p, 1);
  extend_uturn_positive(p, 1);
  extend_uturn_negative(p, 1);
  extend_uturn_positive(p, 1);
  extend_uturn_negative(p, 1);
  extend_uturn_positive(p, 1);
  extend_straight_negative(p, 1);
  extend_uturn_positive(p, 1);
  extend_straight_negative(p, 2);
}

template<typename Path>
void generate_g2_torus(Path& p, std::size_t i)
{
  assert(i<3);
  switch(i)
  {
    case 0: generate_g2_v0_torus(p); break;
    case 1: generate_g2_v1_torus(p); break;
    case 2: generate_g2_v2_torus(p); break;
    default: assert(false);
  }
}

template<typename Path>
void generate_g1_v0_double_torus(Path& p)
{ // 1st generator
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[16]));
  extend_straight_negative(p, 2);
}

template<typename Path>
void generate_g1_v1_double_torus(Path& p)
{ // 1st generator
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[4]));
  extend_straight_negative(p, 2);
}

template<typename Path>
void generate_g1_v2_double_torus(Path& p)
{ // 1st generator
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[2]));
  extend_straight_negative(p, 2);
  p.reverse();
}

template<typename Path>
void generate_g2_v1_double_torus(Path& p)
{ // 2nd generator
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[10]));
  extend_straight_negative(p, 2);
  extend_uturn_negative(p, 1);
}

template<typename Path>
void generate_g1_double_torus(Path& p, std::size_t i)
{
  assert(i<3);
  switch(i)
  {
    case 0: generate_g1_v0_double_torus(p); break;
    case 1: generate_g1_v1_double_torus(p); break;
    case 2: generate_g1_v2_double_torus(p); break;
    default: assert(false);
  }
}

} // namespace CGAL

#endif // CGAL_CREATION_OF_TEST_CASES_FOR_PATHS_H //
// EOF //
