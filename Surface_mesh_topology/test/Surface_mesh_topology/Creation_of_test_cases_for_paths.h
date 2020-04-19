// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_CREATION_OF_TEST_CASES_FOR_PATHS_H
#define CGAL_CREATION_OF_TEST_CASES_FOR_PATHS_H 1

#include<CGAL/Surface_mesh_topology/internal/Path_generators.h>

namespace CGAL {

template<typename Path>
void generate_one_positive_spur(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5]));
  p.extend_straight_positive(6);
  p.extend_positive_turn();
  p.extend_positive_turn(0);
  p.extend_straight_positive(4);
}

template<typename Path>
void generate_one_negative_spur(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5]));
  p.extend_straight_negative(6);
  p.extend_negative_turn();
  p.extend_positive_turn(0);
  p.extend_straight_negative(4);
}

template<typename Path>
void generate_cyclic_spur(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5]));
  p.extend_positive_turn(0);
}

template<typename Path>
void generate_one_positive_bracket(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5]));
  p.extend_straight_positive(3);
  p.extend_positive_turn(3);
  p.extend_positive_turn();
  p.extend_straight_positive(6);
  p.extend_positive_turn();
  p.extend_positive_turn(3);
  p.extend_straight_positive(2);
}

template<typename Path>
void generate_one_negative_bracket(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5]));
  p.extend_straight_negative(3);
  p.extend_negative_turn();
  p.extend_straight_negative(6);
  p.extend_negative_turn();
  p.extend_straight_negative(2);
}

template<typename Path>
void generate_bracket_special1(Path& p, bool reverse)
{ // Case (x, 1, 2^r, 1)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[79]));
  p.extend_positive_turn(1);
  p.extend_straight_positive(7);
  p.extend_positive_turn();

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
  p.extend_positive_turn(1);
  p.extend_straight_positive(10);

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
  p.extend_straight_negative(2);
  p.extend_negative_turn(3);
  p.extend_straight_negative(8);
  p.extend_negative_turn();
  p.extend_straight_negative(5);
  p.extend_negative_turn(3);
  p.extend_straight_negative(3);
}

template<typename Path>
void generate_l_shape_case2(Path& p)
{ // (... x -1 -2^t y ...): here (-2^2 -3 -1 -2^5 -3 -2^3)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[0]));
  p.extend_straight_negative(2);
  p.extend_negative_turn(3);
  p.extend_negative_turn();
  p.extend_straight_negative(5);
  p.extend_negative_turn(3);
  p.extend_straight_negative(3);
}

template<typename Path>
void generate_l_shape_case3(Path& p)
{ // (... x -2^s -1 y ...): here (-2^2 -3 -2^5 -1 -3 -2^3)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[0]));
  p.extend_straight_negative(2);
  p.extend_negative_turn(3);
  p.extend_straight_negative(5);
  p.extend_negative_turn(1);
  p.extend_negative_turn(3);
  p.extend_straight_negative(3);
}

template<typename Path>
void generate_l_shape_case4(Path& p)
{ // (x -2^s -1 -2^t): here (-2^7 -1 -2^3 -4)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[35]));
  p.extend_negative_turn(4);
  p.extend_straight_negative(7);
  p.extend_negative_turn(1);
  p.extend_straight_negative(2);
}

template<typename Path>
void generate_l_shape_case5(Path& p)
{ // (x -1 -2^t): here (-4 -1 -2^12)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[27]));
  p.extend_negative_turn(1);
  p.extend_straight_negative(12);
}

template<typename Path>
void generate_l_shape_case6(Path& p)
{ // (x -2^t -1): here (-4 -2^12 -1)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[13]));
  p.extend_straight_negative(12);
  p.extend_negative_turn(1);
}

template<typename Path>
void generate_l_shape_case7(Path& p)
{ // (-3 -2^s -1 -2^t): here (-2^7 -1 -2^3 -3)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[194]));
  p.extend_straight_negative(7);
  p.extend_negative_turn(1);
  p.extend_straight_negative(3);
}

template<typename Path>
void generate_l_shape_case8(Path& p)
{ // (-2^l): here (-2^20)
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[4]));
  p.extend_straight_negative(19);
}

template<typename Path>
void generate_g1_v0_torus(Path& p)
{ // 1st generator
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[0]));
  p.extend_straight_negative(4);
}

template<typename Path>
void generate_g1_v1_torus(Path& p)
{ // 1st generator v2
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[20]));
  p.extend_straight_negative(4);
}

template<typename Path>
void generate_g1_v2_torus(Path& p)
{ // 1st generator v3
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[20]));
  p.extend_negative_turn(1);
  p.extend_straight_negative(1);
  p.extend_positive_turn(1);
  p.extend_straight_negative(2);
  p.extend_positive_turn(1);
  p.extend_straight_negative(2);
  p.extend_negative_turn(1);
  p.extend_negative_turn(1);
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
  p.extend_positive_turn(1);
  p.extend_positive_turn(1);
  p.extend_positive_turn(1);
}

template<typename Path>
void generate_null_cycle_v1_torus(Path& p)
{ // Empty cycle v2
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[20]));
  p.extend_positive_turn(1);
  p.extend_straight_positive(1);
  p.extend_positive_turn(1);
  p.extend_positive_turn(1);
  p.extend_straight_positive(1);
}

template<typename Path>
void generate_null_cycle_v2_torus(Path& p)
{ // Empty cycle v3
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[20]));
  p.extend_negative_turn(1);
  p.extend_straight_negative(3);
  p.extend_negative_turn(1);
  p.extend_straight_negative(2);
  p.extend_negative_turn(1);
  p.extend_straight_negative(3);
  p.extend_negative_turn(1);
  p.extend_straight_negative(1);
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
  p.extend_straight_positive(4);
}

template<typename Path>
void generate_g2_v1_torus(Path& p)
{ // 2nd generator v2
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[13]));
  p.extend_straight_negative(4);
}

template<typename Path>
void generate_g2_v2_torus(Path& p)
{ // 2nd generator v3
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[0]));
  p.extend_negative_turn(1);
  p.extend_positive_turn(1);
  p.extend_negative_turn(1);
  p.extend_positive_turn(1);
  p.extend_negative_turn(1);
  p.extend_positive_turn(1);
  p.extend_straight_negative(1);
  p.extend_positive_turn(1);
  p.extend_straight_negative(2);
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
  p.extend_straight_negative(2);
}

template<typename Path>
void generate_g1_v1_double_torus(Path& p)
{ // 1st generator
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[4]));
  p.extend_straight_negative(2);
}

template<typename Path>
void generate_g1_v2_double_torus(Path& p)
{ // 1st generator
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[2]));
  p.extend_straight_negative(2);
  p.reverse();
}

template<typename Path>
void generate_g2_v1_double_torus(Path& p)
{ // 2nd generator
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[10]));
  p.extend_straight_negative(2);
  p.extend_negative_turn(1);
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

template<typename Path>
void generate_path(Path& p,
                   const std::initializer_list<int>& v1,
                   const std::initializer_list<int>& v2)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[0]));

  std::vector<int> before(v1);
  for (auto it=before.begin(), itend=before.end(); it!=itend; )
  {
    if (*it<0) { p.extend_negative_turn(-(*it), false); }
    else { p.extend_positive_turn(*it, false); }
    ++it;

    p.extend_straight_positive(*it, false);
    ++it;
  }

  Surface_mesh_topology::internal::create_braket_positive(p, 6, false);

  std::vector<int> after(v2);
  for (auto it=after.begin(), itend=after.end(); it!=itend; )
  {
    if (*it<0) { p.extend_negative_turn(-(*it), false); }
    else { p.extend_positive_turn(*it, false); }
    ++it;

    p.extend_straight_positive(*it, false);
    ++it;
  }

  p.update_is_closed();
}


} // namespace CGAL

#endif // CGAL_CREATION_OF_TEST_CASES_FOR_PATHS_H //
// EOF //
