// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
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
//******************************************************************************

#include "test_meshing_utilities.h"
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/IO/File_tetgen.h>
#include <CGAL/IO/File_binary_mesh_3.h>
#include <CGAL/use.h>

#include <fstream>

template <typename K, typename Concurrency_tag = CGAL::Sequential_tag>
struct Polyhedron_with_features_tester : public Tester<K>
{
  void operator()() const
  {
    typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Gt;
    typedef typename CGAL::Mesh_polyhedron_3<Gt, short>::type Polyhedron;
    typedef CGAL::Polyhedral_mesh_domain_with_features_3<Gt,
                                                         Polyhedron,
                                                         CGAL::Default,
                                                         short> Mesh_domain;

    typedef typename CGAL::Mesh_triangulation_3<
      Mesh_domain,
      typename CGAL::Kernel_traits<Mesh_domain>::Kernel,
      Concurrency_tag>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3 <
      Tr,
      typename Mesh_domain::Corner_index,
      typename Mesh_domain::Curve_index > C3t3;

    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    typedef typename Mesh_criteria::Edge_criteria Edge_criteria;
    typedef typename Mesh_criteria::Facet_criteria Facet_criteria;
    typedef typename Mesh_criteria::Cell_criteria Cell_criteria;

    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    std::cout << "\tSeed is\t"
      << CGAL::get_default_random().get_seed() << std::endl;
    std::ifstream input("data/cube.off");
    Polyhedron polyhedron;
    input >> polyhedron;
    Mesh_domain domain(polyhedron, &CGAL::get_default_random());
    domain.detect_features();

    // non-documented, provided to the FEniCS project
    const std::vector<Polyhedron>& polyhedra = domain.polyhedra();
    CGAL_USE(polyhedra);

    // Set mesh criteria
#ifndef CGAL_MESH_3_VERBOSE
    Edge_criteria edge_criteria(0.2);
    Facet_criteria facet_criteria(30, 0.2, 0.02);
    Cell_criteria cell_criteria(3, 0.2);
#else // a different set of criteria, for the test of CGAL_MESH_3_VERBOSE
    Edge_criteria edge_criteria(0.3);
    Facet_criteria facet_criteria(30, 0.3, 0.03);
    Cell_criteria cell_criteria(3, 0.4);
#endif
    Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        CGAL::parameters::manifold(),
                                        CGAL::parameters::no_exude(),
                                        CGAL::parameters::no_perturb());

    CGAL::remove_far_points_in_mesh_3(c3t3);

    // Verify
    this->verify(c3t3,domain,criteria,
                 Polyhedral_tag()); //, 1099, 1099, 1158, 1158, 4902, 4902);

    std::ofstream out_medit("test-medit.mesh");
    CGAL::IO::output_to_medit(out_medit, c3t3);
    CGAL::IO::output_to_tetgen("test-tetgen", c3t3);
    std::ofstream out_binary("test-binary.mesh.cgal",
                             std::ios_base::out|std::ios_base::binary);
    CGAL::IO::save_binary_file(out_binary, c3t3);
    out_binary.close();
    C3t3 c3t3_bis;
    std::ifstream in_binary("test-binary.mesh.cgal",
                             std::ios_base::in|std::ios_base::binary);
    CGAL::IO::load_binary_file(in_binary, c3t3_bis);
    assert(c3t3_bis.triangulation() == c3t3.triangulation());

  }
};

int main()
{
  Polyhedron_with_features_tester<K_e_i> test_epic;
  std::cerr << "Mesh generation from a polyhedron with edges:\n";
  test_epic();

#ifdef CGAL_LINKED_WITH_TBB
  Polyhedron_with_features_tester<K_e_i, CGAL::Parallel_tag> test_epic_p;
  std::cerr << "Parallel mesh generation from a polyhedron with edges:\n";
  test_epic_p();
#endif

  return EXIT_SUCCESS;
}
