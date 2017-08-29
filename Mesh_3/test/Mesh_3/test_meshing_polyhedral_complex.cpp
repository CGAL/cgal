// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description :
//******************************************************************************

#include "test_meshing_utilities.h"

#include <CGAL/Polyhedral_complex_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/IO/File_binary_mesh_3.h>

#include <fstream>

const char* const filenames[] = {
  "data/patches/patch-01.off",
  "data/patches/patch-13.off",
  "data/patches/patch-20.off",
  "data/patches/patch-21.off",
  "data/patches/patch-23.off",
  "data/patches/patch-30.off",
};

const std::pair<int, int> incident_subdomains[] = {
  std::make_pair(0, 1),
  std::make_pair(1, 3),
  std::make_pair(2, 0),
  std::make_pair(2, 1),
  std::make_pair(2, 3),
  std::make_pair(3, 0),
};


template <typename K, typename Concurrency_tag = CGAL::Sequential_tag>
struct Polyhedral_complex_tester : public Tester<K>
{
  void operator()() const
  {
    typedef typename CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
    typedef CGAL::Polyhedral_complex_mesh_domain_3<K> Mesh_domain;
    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<
      Tr,
      typename Mesh_domain::Corner_index,
      typename Mesh_domain::Curve_segment_index> C3t3;
    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

    //Input
    const std::size_t nb_patches = sizeof(filenames) / sizeof(const char*);
      assert(sizeof(incident_subdomains) ==
        nb_patches * sizeof(std::pair<int, int>));

    std::vector<Polyhedron> patches(nb_patches);
    for (std::size_t i = 0; i < nb_patches; ++i) {
      std::ifstream input(filenames[i]);
      if (!(input >> patches[i])) {
        std::cerr << "Error reading " << filenames[i] << " as a polyhedron!\n";
        return;
      }
    }

    // Create domain
    std::cout << "\tSeed is\t"
      << CGAL::get_default_random().get_seed() << std::endl;
    Mesh_domain domain(patches.begin(), patches.end(),
      incident_subdomains, incident_subdomains + nb_patches,
      &CGAL::get_default_random());

    domain.detect_features(); //includes detection of borders

    //check compilation
    std::set<typename K::Point_3> duplicates;
    domain.merge_duplicated_points(duplicates);
    //end check compilation

    // Mesh criteria
    using namespace CGAL::parameters;
    Mesh_criteria criteria(edge_size = 8,
      facet_angle = 25, facet_size = 8, facet_distance = 0.2,
      cell_radius_edge_ratio = 3, cell_size = 10);

    // Mesh generation
    C3t3 c3t3;

    CGAL::internal::Mesh_3::init_c3t3_with_features(c3t3, domain, criteria,
      true /*nonlinear_growth_of_balls*/);
    domain.add_vertices_to_c3t3_on_patch_without_feature_edges(c3t3);

    CGAL::refine_mesh_3<C3t3>(c3t3, domain, criteria);

    CGAL::remove_far_points_in_mesh_3(c3t3);

    // Verify
    // does not work because of Hausdorff distance
    // which does not mean anything here
    //this->verify(c3t3, domain, criteria, Polyhedral_tag());
  }
};

int main()
{
  Polyhedral_complex_tester<K_e_i> test_epic;
  std::cerr << "Mesh generation from a polyhedral complex:\n";
  test_epic();
  
#ifdef CGAL_LINKED_WITH_TBB
  Polyhedral_complex_tester<K_e_i, CGAL::Parallel_tag> test_epic_p;
  std::cerr << "Parallel mesh generation from a polyhedral complex:\n";
  test_epic_p();
#endif

  return EXIT_SUCCESS;
}
