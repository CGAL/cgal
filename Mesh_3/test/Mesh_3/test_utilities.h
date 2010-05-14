// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef TEST_UTILITIES_H_
#define TEST_UTILITIES_H_


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K_e_i;
typedef CGAL::Exact_predicates_exact_constructions_kernel K_e_e;

namespace CGAL {
namespace details {

template<>
struct Mesh_geom_traits_generator<K_e_e>
{
private:
  typedef Regular_triangulation_euclidean_traits_3<K_e_e> Geom_traits;

public:
  typedef Geom_traits type;
  typedef type Type;
};  // end struct Mesh_geom_traits_generator<...>

} // end namespace details
} // end namespace CGAL


#endif // TEST_UTILITIES_H_
