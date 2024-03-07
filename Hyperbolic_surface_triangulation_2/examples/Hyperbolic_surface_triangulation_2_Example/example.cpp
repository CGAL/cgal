// Copyright (c) 2024
// INRIA Nancy (France), and Université Gustave Eiffel Marne-la-Vallee (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Loïc Dubois, Monique Teillaud

#include <CGAL/Hyperbolic_fundamental_domain_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Hyperbolic_surface_triangulation_2.h>
#include <CGAL/Hyperbolic_surfaces_traits_2.h>
#include <CGAL/Gmpq.h>

#include <time.h>

using namespace CGAL;

typedef Hyperbolic_surfaces_traits_2<Gmpq>                        Traits;
typedef Hyperbolic_fundamental_domain_2<Traits>                   Domain;
typedef Hyperbolic_surface_triangulation_2<Traits>                Triangulation;
typedef Hyperbolic_fundamental_domain_factory_2<Traits>           Factory;

int main(){
  Factory factory = Factory(time(NULL));
  Domain domain0, domain;
  Triangulation triangulation;
  std::string filename = "./data/domain";

    // Generate the domain
    std::cout << "generating the domain " << std::endl;
    domain0 = factory.generate_domain_g2();

    // Save the domain
    std::cout << "saving the domain " << std::endl;
    std::ofstream output_file (filename.c_str());
    output_file << domain0;
    output_file.close();

    // Load the domain
    std::cout << "loading the domain" << std::endl;
    std::ifstream input_file (filename.c_str());
    input_file >> domain;
    input_file.close();

    // Triangulate the domain
    std::cout << "triangulating the domain " << std::endl;
    triangulation = Triangulation(domain);

    // Save the resulting triangulation
    std::cout << "saving the input triangulation " << std::endl;
    output_file = std::ofstream ("./data/input triangulation");
    output_file << triangulation;
    output_file.close();

    // Delaunay flip the triangulation
    std::cout << "Delaunay flipping the triangulation " << std::endl;
    triangulation.make_delaunay();

    // Save the resulting triangulation
    std::cout << "saving the output triangulation " << std::endl;
    output_file = std::ofstream ("./data/output triangulation");
    output_file << triangulation;
    output_file.close();

    // Load the resulting triangulation
    std::cout << "loading the resulting triangulation" << std::endl;
    Triangulation triangulation2;
    input_file = std::ifstream ("./data/output triangulation");
    input_file >> triangulation2;
    input_file.close();

    std::cout << triangulation2 << std::endl;

  return 0;
}
