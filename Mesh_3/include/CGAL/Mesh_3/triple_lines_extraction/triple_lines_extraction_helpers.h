// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot, Jane Tournois
//
//******************************************************************************
//
//******************************************************************************


#ifndef CGAL_MESH_3_TRIPLE_LINES_HELPERS_H
#define CGAL_MESH_3_TRIPLE_LINES_HELPERS_H

#include <CGAL/license/Mesh_3.h>

namespace CGAL
{
namespace Mesh_3
{
namespace internal
{

  template<typename Word_type, bool b = (sizeof(Word_type) > 1)>
  struct Color_transformation_helper
  {
    using type = std::unordered_map<Word_type, std::uint8_t>;
    static void reset(type& t)
    {
      t.clear();
    }
    static bool is_valid(const type& t, const Word_type& w)
    {
      return t.find(w) != t.end();
    }
  };
  template<typename Word_type>
  struct Color_transformation_helper<Word_type, false>
  {
    using type = std::array<std::uint8_t, (1 << (sizeof(Word_type) * 8))>;
    static void reset(type& t)
    {
      std::fill(t.begin(), t.end(), 8/*invalid_word*/);
    }
    static bool is_valid(const type& t, const Word_type& w)
    {
      return t[w] != 8;/*invalid_word*/
    }
  };

  template<typename T, int N>
  void debug_cerr(const char* title, const std::array<T, N>& tab)
  {
    std::cerr << title << " [";
    for (const T& t : tab)
      std::cout << (int)(t) << " ";
    std::cout << "]" << std::endl;
  }

}//end namespace internal
}//end namespace Mesh_3
}//end namespace CGAL

#endif // CGAL_MESH_3_TRIPLE_LINES_HELPERS_H
