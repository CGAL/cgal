// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
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
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

template <typename K>
struct Polyhedron_with_features_tester : public Tester<K>
{
  void operator()() const
  {
    typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Gt;
    typedef CGAL::Polyhedral_mesh_domain_with_features_3<Gt> Mesh_domain;
    
    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3 <
      Tr,
      typename Mesh_domain::Corner_index,
      typename Mesh_domain::Curve_segment_index > C3t3;
    
    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    typedef typename Mesh_criteria::Edge_criteria Edge_criteria;
    typedef typename Mesh_criteria::Facet_criteria Facet_criteria;
    typedef typename Mesh_criteria::Cell_criteria Cell_criteria;
    
    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    Mesh_domain domain("data/cube.off");
    domain.detect_features();
    
    // Set mesh criteria
    Edge_criteria edge_criteria(0.2);
    Facet_criteria facet_criteria(30, 0.2, 0.02);
    Cell_criteria cell_criteria(3, 0.2);
    Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        CGAL::parameters::no_exude(),
                                        CGAL::parameters::no_perturb());
    
    // Verify
    this->verify(c3t3,domain,criteria); //, 1099, 1099, 1158, 1158, 4902, 4902);
  }
};

int main()
{
  Polyhedron_with_features_tester<K_e_i> test_epic;
  std::cerr << "Mesh generation from a polyhedron with edges:\n";
  test_epic();

  return EXIT_SUCCESS;
}
