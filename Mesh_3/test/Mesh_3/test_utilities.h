// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_MESH_3_TEST_TEST_UTILITIES_H
#define CGAL_MESH_3_TEST_TEST_UTILITIES_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/disable_warnings.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K_e_i;
typedef CGAL::Exact_predicates_exact_constructions_kernel K_e_e;

namespace CGAL {
namespace details {

template<>
struct Mesh_geom_traits_generator<K_e_e>
{
private:
  typedef K_e_e Geom_traits;

public:
  typedef Geom_traits type;
  typedef type Type;
};  // end struct Mesh_geom_traits_generator<...>

} // end namespace details
} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_3_TEST_TEST_UTILITIES_H
