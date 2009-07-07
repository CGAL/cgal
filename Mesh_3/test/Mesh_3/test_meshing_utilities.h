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
// File Description : Test meshing utilities.
//******************************************************************************

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>

#include "test_utilities.h"

#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include <vector>

// IO
#include <fstream>
#include <iostream>


template <typename K>
struct Tester
{
  template<typename C3t3, typename Domain, typename Criteria>
  void verify(C3t3& c3t3,
              const Domain& domain,
              const Criteria& criteria,
              const unsigned int min_vertices_expected = 0,
              const unsigned int max_vertices_expected = (unsigned int)(-1),
              const unsigned int min_facets_expected = 0,
              const unsigned int max_facets_expected = (unsigned int)(-1),
              const unsigned int min_cells_expected = 0,
              const unsigned int max_cells_expected = (unsigned int)(-1) ) const
  {
    typedef typename C3t3::size_type size_type;

    // Store mesh properties
    size_type v = c3t3.triangulation().number_of_vertices();
    size_type f = c3t3.number_of_facets();
    size_type c = c3t3.number_of_cells();

    // Verify
    verify_c3t3(c3t3,
                min_vertices_expected,
                max_vertices_expected,
                min_facets_expected,
                max_facets_expected,
                min_cells_expected,
                max_cells_expected);

    // Refine again and verify nothing changed
    std::cerr << "\tRefining again...\n";
    refine_mesh_3(c3t3,domain,criteria,false);
    verify_c3t3(c3t3,v,v,f,f,c,c);
  }

  template<typename C3t3>
  void verify_c3t3(const C3t3& c3t3,
                   const unsigned int min_vertices_expected = 0,
                   const unsigned int max_vertices_expected = (unsigned int)(-1),
                   const unsigned int min_facets_expected = 0,
                   const unsigned int max_facets_expected = (unsigned int)(-1),
                   const unsigned int min_cells_expected = 0,
                   const unsigned int max_cells_expected = (unsigned int)(-1)) const
  {
    //-------------------------------------------------------
    // Verifications
    //-------------------------------------------------------
    std::cerr << "\t\tNumber of cells: " << c3t3.number_of_cells() << "\n";
    std::cerr << "\t\tNumber of facets: " << c3t3.number_of_facets() << "\n";
    std::cerr << "\t\tNumber of vertices: " << c3t3.number_of_vertices() << "\n";

    assert(min_vertices_expected <= c3t3.number_of_vertices());
    assert(max_vertices_expected >= c3t3.number_of_vertices());

    assert(min_facets_expected <= c3t3.number_of_facets());
    assert(max_facets_expected >= c3t3.number_of_facets());

    assert(min_cells_expected <= c3t3.number_of_cells());
    assert(max_cells_expected >= c3t3.number_of_cells());
  }
};

