// Copyright (c) 2025 GeometryFactory, France.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description :
//******************************************************************************

#include <CGAL/SMDS_3/io_signature.h>
#include "test_meshing_utilities.h"
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/SMDS_3/Dump_c3t3.h>

#include <CGAL/disable_warnings.h>

#include <type_traits>
#include <string>

template <typename K>
struct Polyhedron_tester : public Tester<K>
{
  void polyhedron(const std::string& filename) const
  {
    using GT = K;
    using Polyhedron = CGAL::Polyhedron_3<GT>;
    using Mesh_domain = CGAL::Polyhedral_mesh_domain_3<Polyhedron, GT>;
    using Tr = typename CGAL::Mesh_triangulation_3<Mesh_domain, GT>::type;
    using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;
    using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

    Polyhedron polyhedron;
    std::ifstream input(filename);
    input >> polyhedron;
    input.close();

    auto rng = &CGAL::get_default_random();
    const auto seed = rng->get_seed();

    std::cout << "\tSeed is\t" << seed << std::endl;
    Mesh_domain domain(polyhedron, rng);

    // Set mesh criteria - only facet criteria here
    Mesh_criteria criteria(CGAL::parameters::facet_angle(30)
                                            .facet_size(0.2)
                                            .facet_distance(0.02));

    // Mesh generation with surface_only() option
    C3t3 c3t3_surface_only = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                                     CGAL::parameters::surface_only());

    assert(c3t3_surface_only.number_of_facets_in_complex() > 0);
    assert(c3t3_surface_only.number_of_cells_in_complex() == 0);

    CGAL::dump_c3t3(c3t3_surface_only, "c3t3_surface_only.binary.cgal");

    std::cout << "\t With surface_only() option:" << std::endl;
    std::cout << "\t\t Nb of facets = " << c3t3_surface_only.number_of_facets_in_complex() << std::endl;
    std::cout << "\t\t Nb of cells = " << c3t3_surface_only.number_of_cells_in_complex() << std::endl;

    // Mesh generation without surface_only() option,
    // i.e. volume meshing, with no cell criteria
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        CGAL::parameters::no_exude().no_perturb());

    assert(c3t3.number_of_facets_in_complex() > 0);
    if(is_closed(polyhedron)) {
      assert(c3t3.number_of_cells_in_complex() > 0);
    } else {
      assert(c3t3.number_of_cells_in_complex() == 0);
    }

    CGAL::dump_c3t3(c3t3, "c3t3_no_surface_only.binary.cgal");

    std::cout << "\t Without surface_only() option:" << std::endl;
    std::cout << "\t\t Nb of facets = " << c3t3.number_of_facets_in_complex() << std::endl;
    std::cout << "\t\t Nb of cells = " << c3t3.number_of_cells_in_complex() << std::endl;
  }
};

int main()
{
  Polyhedron_tester<K_e_i> test_epic;

  std::cerr << "Remeshing from a closed polyhedron (sphere.off):\n";
  const std::string sphere = CGAL::data_file_path("meshes/sphere.off");
  test_epic.polyhedron(sphere);

  std::cerr << "Remeshing from a non-closed polyhedron (lion.off):\n";
  const std::string lion = CGAL::data_file_path("meshes/lion.off");
  test_epic.polyhedron(lion);

  return EXIT_SUCCESS;
}
