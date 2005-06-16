// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source: 
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Steve OUDOT


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Complex_2_in_triangulation_vertex_base_3.h>
#include <CGAL/Complex_2_in_triangulation_surface_mesh_cell_base_3.h>
#include <CGAL/Complex_2_in_triangulation_vertex_base_3.h>
#include <CGAL/Robust_circumcenter_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesher/Surface_mesher.h>
#include <CGAL/Surface_mesher/Surface_mesher_regular_edges.h>
#include <CGAL/Surface_mesher/Surface_mesher_regular_edges_without_boundary.h>
#include <CGAL/Surface_mesher/Surface_mesher_manifold.h>
#include <CGAL/Surface_mesher/Criteria/Standard_criteria.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

// #include <CGAL/Surface_mesher/Oracles/Implicit_oracle.h>
// #include "implicit_function.h"
#include <CGAL/Surface_mesher/Oracles/Polyhedral.h>

#include <fstream>




/////////////// Types /////////////// 

struct K2 : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Robust_circumcenter_traits_3<K2>  K;
typedef CGAL::Triangulation_vertex_base_3<K> Vb;
typedef CGAL::Complex_2_in_triangulation_vertex_base_3<K, Vb> Vb2;
typedef CGAL::Complex_2_in_triangulation_surface_mesh_cell_base_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb2, Cb> Tds;

typedef CGAL::Delaunay_triangulation_3<K, Tds> Del;
typedef CGAL::Complex_2_in_triangulation_3_surface_mesh<Del> C2t3;
// typedef Function <K::FT> Func;
// typedef CGAL::Surface_mesher::Implicit_oracle<K, Func> Oracle;
typedef CGAL::Surface_mesher::Polyhedral <Del> Oracle;
typedef CGAL::Surface_mesher::Refine_criterion<Del> Criterion;
typedef CGAL::Surface_mesher::Standard_criteria <Criterion > Criteria;

typedef CGAL::Surface_mesher::Surface_mesher<Del, Oracle, Criteria> SM;
typedef CGAL::Surface_mesher::Surface_mesher_regular_edges<Del, Oracle, Criteria> SMRE;
typedef CGAL::Surface_mesher::Surface_mesher_regular_edges_without_boundary<Del, Oracle, Criteria> SMREWB;
typedef CGAL::Surface_mesher::Surface_mesher_manifold<Del, Oracle, Criteria> SMM;
typedef CGAL::Surface_mesher::Surface_mesher_regular_edges_without_boundary_base<Del, Oracle, Criteria> SMREWBB;
typedef CGAL::Surface_mesher::Surface_mesher_manifold<Del, Oracle, Criteria, 
CGAL::Surface_mesher::Surface_mesher_manifold_base <Del, Oracle, Criteria, SMREWBB> 
> SMMWB;


typedef SM Surface;  // basic mesher
//typedef SMRE Surface;
//typedef SMREWB Surface;
// typedef SMM Surface;  // manifold with boundary
//typedef SMMWB Surface;  // manifold without boundary



/////////////// Main function /////////////// 

int main(int argc, char **argv) {

  if (argc < 2) {
    std::cout << "Usage : " << argv[0] << " <off_file>" << std::endl;
    exit(0);
  }


  // Function
//   Func F;

  // Oracle (NB: parity oracle is toggled)
//   Oracle O (F, K::Point_3 (0,0,0), 4, 1e-6, true);  
  std::ifstream is (argv[1]);
  Oracle O (is);
  is.close ();

  // 3D-Delaunay triangulation
  Del T;


  // Initial point sample
  //Oracle::Points initial_point_sample = O.random_points (10);
  //  Oracle::Points initial_point_sample = O.random_points (20);
  Oracle::Points initial_point_sample = O.random_points (50);
  typedef Del::Point Point;
  T.insert (initial_point_sample.begin(), initial_point_sample.end());
  
  // Meshing criteria
  CGAL::Surface_mesher::Curvature_size_criterion<Del> c_s_crit (0.01);
  CGAL::Surface_mesher::Uniform_size_criterion<Del> u_s_crit (10000);
  CGAL::Surface_mesher::Aspect_ratio_criterion<Del> a_r_crit (30);
  std::vector<Criterion*> crit_vect;
  crit_vect.push_back (&a_r_crit);
  crit_vect.push_back (&u_s_crit);
  crit_vect.push_back (&c_s_crit);
  //std::vector<Criterion*> crit_vect(1);
  //crit_vect[0] = &u_s_crit;
  //crit_vect[1] = &a_r_crit;
  Criteria C (crit_vect);

  std::cout << "Initial number of points: " << T.number_of_vertices() 
            << std::endl;
  
  // 2D-complex in 3D-Delaunay triangulation
  C2t3 Co2 (T);
  // Surface meshing
  Surface mesher (T, Co2, O, C);
//   Surface mesher (T, O, a_r_crit);
  mesher.refine_mesh(true);
  
  std::cout << "Final number of points: " << T.number_of_vertices() 
            << std::endl;

  // Output
  if (argc >= 3) {
    std::cout << "Writing output to file " << argv[2] << "...";
    std::ofstream os(argv[2]);
    output_surface_facets_to_off (os, T);
    os.close();
    std::cout << " done\n";
  }
}
