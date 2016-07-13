// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#include <CGAL/Mesh_3/io_signature.h>
#include "test_meshing_utilities.h"
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <boost/type_traits/is_same.hpp>

#include <CGAL/Mesh_3/Dump_c3t3.h>

template <typename K, typename Concurrency_tag = CGAL::Sequential_tag>
struct Polyhedron_tester : public Tester<K>
{
  void polyhedron() const
  {
    typedef K Gt;
    typedef CGAL::Polyhedron_3<Gt> Polyhedron;
    typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Gt> Mesh_domain;

    CGAL_static_assertion((boost::is_same<
                            typename Mesh_domain::Surface_patch_index,
                            std::pair<int, int>
                           >::value));

    typedef typename CGAL::Mesh_triangulation_3<
      Mesh_domain,
      typename CGAL::Kernel_traits<Mesh_domain>::Kernel,
      Concurrency_tag>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    typedef typename Mesh_criteria::Facet_criteria Facet_criteria;
    typedef typename Mesh_criteria::Cell_criteria Cell_criteria;

    typedef typename Mesh_domain::Surface_patch_index Surface_patch_index;

    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    Polyhedron polyhedron;
    std::ifstream input("data/sphere.off");
    input >> polyhedron;
    input.close();

    std::cout << "\tSeed is\t"
      << CGAL::get_default_random().get_seed() << std::endl;
    Mesh_domain domain(polyhedron, &CGAL::get_default_random());

    // Set mesh criteria
    Facet_criteria facet_criteria(30, 0.2, 0.02);
    Cell_criteria cell_criteria(2, 0.2);
    Mesh_criteria criteria(facet_criteria, cell_criteria);

    // Mesh generation
    C3t3 c3t3;
    typename Polyhedron::Point_iterator end = polyhedron.points_begin();
    int i=0;
    while ( i++ < 4 ) { ++end; }

    c3t3.insert_surface_points(polyhedron.points_begin(),
                               end,
                               domain.index_from_surface_patch_index(Surface_patch_index(0,1)));

    CGAL::refine_mesh_3(c3t3, domain, criteria,
                        CGAL::parameters::no_exude(),
                        CGAL::parameters::no_perturb());

    CGAL::remove_far_points_in_mesh_3(c3t3);
    
    // Verify
    double vol = 0.479171765761454;
    this->verify_c3t3_volume(c3t3, vol*0.95, vol*1.05);
#ifdef CGAL_LINKED_WITH_TBB
    // Parallel
    if (boost::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value)
    {
      this->verify(c3t3, domain, criteria, Polyhedral_tag(),
                   110, 140, 190, 235, 350, 420); 
    }
    else
#endif //CGAL_LINKED_WITH_TBB
    {
      this->verify(c3t3, domain, criteria, Polyhedral_tag(),
                   119, 121, 200, 204, 350, 360);  
    }

    // test the dump function
    CGAL::dump_c3t3(c3t3, "test_meshing_polyhedron-out");
  }
};

int main()
{
  Polyhedron_tester<K_e_i> test_epic;
  std::cerr << "Mesh generation from a polyhedron:\n";
  test_epic.polyhedron();

#ifdef CGAL_LINKED_WITH_TBB
  Polyhedron_tester<K_e_i, CGAL::Parallel_tag> test_epic_parallel;
  std::cerr << "Mesh generation from a polyhedron using Parallel_tag:\n";
  test_epic_parallel.polyhedron();
#else
  std::cerr << "TBB is not installed, or not configured."
            << "The parallel version cannot be tested.\n";
#endif

  return EXIT_SUCCESS;
}
