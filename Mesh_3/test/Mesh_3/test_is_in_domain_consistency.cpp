// Copyright (c) 2026 GeometryFactory (France).
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

#include "test_meshing_utilities.h"
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/disable_warnings.h>

#include <type_traits>


int main()
{
  typedef K_e_i GT; //Epick
  typedef CGAL::Polyhedron_3<GT> Polyhedron;
  typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, GT> Mesh_domain;

  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

  Polyhedron polyhedron;
  std::ifstream input("data/cube-minus-sphere.off");
  input >> polyhedron;
  input.close();

  std::cout << "\tSeed is\t" << CGAL::get_default_random().get_seed() << std::endl;
  Mesh_domain domain(polyhedron, &CGAL::get_default_random());

  Mesh_criteria criteria(CGAL::parameters::facet_distance(0.01).cell_size(0.1));

  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);//, CGAL::parameters::no_exude().no_perturb());

  auto in_domain = domain.is_in_domain_object();

  const auto& tr = c3t3.triangulation();
  for (auto c : tr.finite_cell_handles())
  {
    const auto subdomain_in_cell = c->subdomain_index();
    const auto cc = tr.dual(c);
    
    std::cout << "\ncc = " << cc << std::endl;
    std::cout << "tet = "
              << c->vertex(0)->point().point() << "\n"
              << c->vertex(1)->point().point() << "\n"
              << c->vertex(2)->point().point() << "\n"
              << c->vertex(3)->point().point() << "\n";

    const auto subdomain_from_oracle = in_domain(tr.dual(c));
    if(subdomain_from_oracle == std::nullopt)
      assert(subdomain_in_cell == 0);
    else
      assert(subdomain_in_cell == subdomain_from_oracle);
  }

  return EXIT_SUCCESS;
}