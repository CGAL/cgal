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
#include <CGAL/Implicit_mesh_domain_3.h>

template <typename K>
struct Implicit_tester : public Tester<K>
{
  typedef typename K::Point_3 Point;
  typedef typename K::FT FT;
  static FT sphere_function (const Point& p)
  {
    const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
    return x2+y2+z2-1;
  }
  
  void implicit() const
  {
    
    typedef FT (Function)(const Point&);
    
    typedef CGAL::Implicit_mesh_domain_3<Function, K> Mesh_domain;
    
    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
    
    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    typedef typename Mesh_criteria::Facet_criteria Facet_criteria;
    typedef typename Mesh_criteria::Cell_criteria Cell_criteria;
    
    typedef typename K::Sphere_3 Sphere_3;
    
    typedef typename Mesh_domain::Surface_index Surface_index;
    
    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    Mesh_domain domain(Implicit_tester<K>::sphere_function,
                       Sphere_3(CGAL::ORIGIN,2.));
    
    // Set mesh criteria
    Facet_criteria facet_criteria(0, 0, 0.3);
    Cell_criteria cell_criteria(0, 0.5);
    Mesh_criteria criteria(facet_criteria, cell_criteria);
    
    std::vector<Point> initial_points;
    initial_points.push_back(Point(1,0,0));
    initial_points.push_back(Point(0,1,0));
    initial_points.push_back(Point(0,0,1));
    initial_points.push_back(Point(-1,0,0));
    initial_points.push_back(Point(0,-1,0));
    initial_points.push_back(Point(0,0,-1));
    
    // Mesh generation
    C3t3 c3t3;
    c3t3.insert_surface_points(initial_points.begin(),
                               initial_points.end(),
                               domain.index_from_surface_index(Surface_index(0,1)));
    
    CGAL::refine_mesh_3(c3t3, domain, criteria,
                        CGAL::parameters::no_exude(),
                        CGAL::parameters::no_perturb());
    
    // Verify
    this->verify(c3t3,domain,criteria,50,58,80,90);
  }
};


int main()
{
  std::cerr << "TESTING WITH Exact_predicates_inexact_constructions_kernel...\n";
  Implicit_tester<K_e_i> test_epic;
  std::cerr << "Mesh generation from an implicit function:\n";
  test_epic.implicit();
  
//  std::cerr << "TESTING WITH Filtered_kernel<Cartesian<float> > kernel...\n";
//  Implicit_tester<Filtered_kernel<CGAL::Cartesian<float> > > test_cf;
//  std::cerr << "Mesh generation from an implicit function:\n";
//  test_cf.implicit();
  
  return EXIT_SUCCESS;

   
  //  std::cerr << "TESTING WITH Filtered_kernel<Simple_cartesian<float> > kernel...\n";
  //  Tester<Filtered_kernel<CGAL::Simple_cartesian<float> > > test_scf;
  //  std::cerr << "Mesh generation from a polyhedron:\n";
  //  test_scf.polyhedron();
  //  std::cerr << "Mesh generation from an implicit function:\n";
  //  test_scf.implicit();
  //
  //  std::cerr << "TESTING WITH Filtered_kernel<Cartesian<float> > kernel...\n";
  //  Tester<Filtered_kernel<CGAL::Cartesian<float> > > test_cf;
  //  std::cerr << "Mesh generation from a polyhedron:\n";
  //  test_cf.polyhedron();
  //  std::cerr << "Mesh generation from an implicit function:\n";
  //  test_cf.implicit();
  //
  //  std::cerr << "TESTING WITH Filtered_kernel<Cartesian<double> > kernel...\n";
  //  Tester<Filtered_kernel<CGAL::Cartesian<double> > > test_cd;
  //  std::cerr << "Mesh generation from a polyhedron:\n";
  //  test_cd.polyhedron();
  //  std::cerr << "Mesh generation from an implicit function:\n";
  //  test_cd.implicit();
  
  //  std::cerr << "\nTESTING WITH Exact_predicates_exact_constructions_kernel...\n";
  //  Tester<K_e_e> test_epec;
  //  std::cerr << "Mesh generation from a polyhedron:\n";
  //  test_epec.polyhedron();
  //  std::cerr << "Mesh generation from an implicit function:\n";
  //  test_epec.implicit();
  //  std::cerr << "Mesh generation from a 3D image:\n";
  //  test_epec.image();
  
}
