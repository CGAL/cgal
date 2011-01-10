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
//******************************************************************************

#include "test_meshing_utilities.h"
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

template <typename K>
struct Polyhedron_tester : public Tester<K>
{
  void polyhedron() const
  {
    typedef K Gt;
    typedef CGAL::Polyhedron_3<Gt> Polyhedron;
    typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Gt> Mesh_domain;
    
    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
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
    
    Mesh_domain domain(polyhedron);
    
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
    
    // Verify
    this->verify(c3t3,domain,criteria, 119, 121, 200, 204, 350, 360);
  }
};

int main()
{
  Polyhedron_tester<K_e_i> test_epic;
  std::cerr << "Mesh generation from a polyhedron:\n";
  test_epic.polyhedron();
  
  return EXIT_SUCCESS;
}
