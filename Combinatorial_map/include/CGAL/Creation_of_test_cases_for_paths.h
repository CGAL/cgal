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
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5])); // 6th dart of the map
  extend_straight_positive(p, 6);
  extend_uturn_positive(p);
  extend_uturn_half_turn(p);
  extend_straight_positive(p, 4);
}

template<typename Path>
void generate_one_negative_spur(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5])); // 6th dart of the map
  extend_straight_negative(p, 6);
  extend_uturn_negative(p);
  extend_uturn_half_turn(p);
  extend_straight_negative(p, 4);
}

template<typename Path>
void generate_cyclic_spur(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5])); // 6th dart of the map
  extend_uturn_half_turn(p);
}

template<typename Path>
void generate_one_positive_bracket(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5])); // 6th dart of the map
  extend_straight_positive(p, 3);
  extend_uturn_negative(p);
  extend_uturn_positive(p);
  extend_straight_positive(p, 6);
  extend_uturn_positive(p);
  extend_uturn_negative(p);
  extend_straight_positive(p, 2);
}

template<typename Path>
void generate_one_negative_bracket(Path& p)
{
  p.clear();
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[5])); // 6th dart of the map
  extend_straight_negative(p, 3);
  extend_uturn_positive(p);
  extend_uturn_negative(p);
  extend_straight_negative(p, 6);
  extend_uturn_negative(p);
  extend_uturn_positive(p);
  extend_straight_negative(p, 2);
}

template<typename Path>
void generate_bracket_special1(Path& p, bool reverse)
{ // Case (x, 1, 2^r, 1)
  p.clear();
  // p.push_back(p.get_map().template beta<2>(p.get_map().darts().iterator_to(p.get_map().darts()[5]))); // 6th dart of the map
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[91])); // 6th dart of the map
  extend_uturn_positive(p, 1);
  extend_straight_positive(p, 8);
  extend_uturn_positive(p);

  if (reverse)
  {
    p.reverse();
    std::cout<<"SPECIAL CASE 1 (x, -1, -2^r, -1): ";
    p.display_negative_turns(); std::cout<<std::endl;
  }
  else
  {
    std::cout<<"SPECIAL CASE 1 (x, 1, 2^r, 1): ";
    p.display_positive_turns(); std::cout<<std::endl;
  }
}

template<typename Path>
void generate_positive_bracket_special1(Path& p)
{ generate_bracket_special1(p, false); }


template<typename Path>
void generate_negative_bracket_special1(Path& p)
{ generate_bracket_special1(p, true); }

} // namespace CGAL

#endif // CGAL_CREATION_OF_TEST_CASES_FOR_PATHS_H //
// EOF //
