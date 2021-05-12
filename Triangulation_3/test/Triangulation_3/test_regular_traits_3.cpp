// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// $Date$
//
//
// Author(s)     : Mariette Yvinec

#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

#include <cassert>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/_test_cls_regular_euclidean_traits_3.h>

//needs exact constructions as well to test the traits class

typedef CGAL::Exact_predicates_exact_constructions_kernel K;

// Explicit instantiation of the whole class :
template class CGAL::Regular_triangulation_euclidean_traits_3<K, K::FT>;

int main()
{
  typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Traits;
  _test_cls_regular_euclidean_traits_3(Traits() );
  std::cerr << "done"<< std::endl;

  return 0;
}
