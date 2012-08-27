// Copyright (c) 2010  GeometryFactory Sarl (France).
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
// Author(s)     : Laurent Rineau <Laurent.Rineau__CGAL@normalesup.org>
// 

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_structural_filtering_traits.h>

#include <cassert>

using namespace CGAL;

struct K : public Epick {};

struct K2 : public Epick {};

namespace CGAL {
template <> struct Triangulation_structural_filtering_traits<K2> {
  typedef Tag_true Use_structural_filtering_tag;
};
}

int main()
{
  assert(Triangulation_structural_filtering_traits<Epeck>::Use_structural_filtering_tag::value == true);
  assert(Triangulation_structural_filtering_traits<Epick>::Use_structural_filtering_tag::value == true);
  assert(Triangulation_structural_filtering_traits<K>::Use_structural_filtering_tag::value == false);
  assert(Triangulation_structural_filtering_traits<K2>::Use_structural_filtering_tag::value == true);
}
