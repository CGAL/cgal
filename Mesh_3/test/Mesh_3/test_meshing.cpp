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
// File Description : Tests meshing.
//
//******************************************************************************


//#define CGAL_MESHER_3_SCAN_VERBOSE
//#define CGAL_MESH_3_DEBUG_CELL_CRITERIA
//#define CGAL_SURFACE_MESHER_DEBUG_CRITERIA
//#define CGAL_MESH_3_VERBOSE
//#define CGAL_MESH_3_DEBUG_FACET_CRITERIA

#include <CGAL/Bbox_3.h>

#include <CGAL/AABB_intersections.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>

#include "test_utilities.h"

#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include <CGAL/Image_3.h>
#include <vector>

// IO
#include <fstream>
#include <iostream>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/File_medit.h>

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>


template <typename K>
struct Tester
{
  void polyhedron() const
  {
    typedef CGAL::Polyhedron_3<K> Polyhedron;
    typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    typedef typename Mesh_criteria::Facet_criteria Facet_criteria;
    typedef typename Mesh_criteria::Cell_criteria Cell_criteria;

    typedef typename Mesh_domain::Surface_index Surface_index;

    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    Polyhedron polyhedron;
    std::ifstream input("data/cube.off");
    input >> polyhedron;
    input.close();

    Mesh_domain domain(polyhedron);

    // Set mesh criteria
    Facet_criteria facet_criteria(25, 0.8, 1);
    Cell_criteria cell_criteria(5, 1.2);
    Mesh_criteria criteria(facet_criteria, cell_criteria);

    // Mesh generation
    C3t3 c3t3;
    c3t3.insert_surface_points(polyhedron.points_begin(),
                               polyhedron.points_end(),
                               domain.index_from_surface_index(Surface_index(0,1)));

    CGAL::refine_mesh_3(c3t3, domain, criteria);

    // Verify
    verify(c3t3,domain,criteria,26,26,48,48);
  }

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
    Mesh_domain domain(Tester<K>::sphere_function, Sphere_3(CGAL::ORIGIN,2.));

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

    CGAL::refine_mesh_3(c3t3, domain, criteria);

    // Verify
    verify(c3t3,domain,criteria,50,58,80,90);
  }

  void image() const
  {
    typedef CGAL::Image_3 Image;
    typedef CGAL::Labeled_image_mesh_domain_3<Image, K> Mesh_domain;

    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    typedef typename Mesh_criteria::Facet_criteria Facet_criteria;
    typedef typename Mesh_criteria::Cell_criteria Cell_criteria;

    typedef typename Mesh_domain::Surface_index Surface_index;

    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    Image image;
    image.read("data/liver_kidney_gallbladder.inr");
    Mesh_domain domain(image,1e-9);

    // Set mesh criteria
    Facet_criteria facet_criteria(25, 20*image.vx(), 5*image.vx());
    Cell_criteria cell_criteria(4, 25*image.vx());
    Mesh_criteria criteria(facet_criteria, cell_criteria);

    // Mesh generation
    C3t3 c3t3 = make_mesh_3<C3t3>(domain,criteria);

    // Verify
    verify(c3t3,domain,criteria);
  }

  template<typename C3t3, typename Domain, typename Criteria>
  void verify(C3t3& c3t3,
              const Domain& domain,
              const Criteria& criteria,
              const unsigned int min_vertices_expected = 0,
              const unsigned int max_vertices_expected = -1,
              const unsigned int min_facets_expected = 0,
              const unsigned int max_facets_expected = -1,
              const unsigned int min_cells_expected = 0,
              const unsigned int max_cells_expected = -1) const
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
    refine_mesh_3(c3t3,domain,criteria);
    verify_c3t3(c3t3,v,v,f,f,c,c);
  }

  template<typename C3t3>
  void verify_c3t3(const C3t3& c3t3,
                   const unsigned int min_vertices_expected = 0,
                   const unsigned int max_vertices_expected = -1,
                   const unsigned int min_facets_expected = 0,
                   const unsigned int max_facets_expected = -1,
                   const unsigned int min_cells_expected = 0,
                   const unsigned int max_cells_expected = -1) const
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






int main()
{
  std::cerr << "TESTING WITH Exact_predicates_inexact_constructions_kernel...\n";
  Tester<K_e_i> test_epic;
  std::cerr << "Mesh generation from a polyhedron:\n";
  test_epic.polyhedron();
  std::cerr << "Mesh generation from an implicit function:\n";
  test_epic.implicit();
  std::cerr << "Mesh generation from a 3D image:\n";
  test_epic.image();

//  std::cerr << "Test trilinear interpolation:\n";
//  test_epic.test();

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

};


// BBox test in AABB_tree...
//FT inv_dir_x(1e9);
//FT inv_dir_y(1e9);
//FT inv_dir_z(1e9);
//
//if ( 0 != direction.x() )
//  inv_dir_x = (FT)1.0/direction.x();
//if ( 0 != direction.y() )
//  inv_dir_y = (FT)1.0/direction.y();
//if ( 0 != direction.z() )
//  inv_dir_y = (FT)1.0/direction.z();
//
//const Vector inv_direction(inv_dir_x, inv_dir_y, inv_dir_z);
