#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Complex_2_in_triangulation_cell_base_3.h>
#include <CGAL/Mesh_3/Complex_2_in_triangulation_cell_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Implicit_surfaces_mesher_3.h>
#include <CGAL/Chew_4_surfaces/Criteria/Standard_criteria.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <CGAL/Chew_4_surfaces/Oracles/Implicit_oracle.h>

#include <CGAL/Mesh_criteria_3.h>

#include "implicit_function.h"

#include <iostream>
#include <fstream>

/////////////// Types /////////////// 

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Triangulation_vertex_base_3<K> Vb;
typedef CGAL::Complex_2_in_triangulation_cell_base_3<K> CCb;
typedef CGAL::Mesh_3::Complex_2_in_triangulation_cell_base_3<K, CCb> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Del;
typedef Function <K::FT> Func;
typedef CGAL::Chew_4_surfaces::Implicit_oracle<K, Func> Oracle;
typedef CGAL::Chew_4_surfaces::Refine_criterion<Del> Criterion;
typedef CGAL::Chew_4_surfaces::Standard_criteria <Criterion > Criteria;
typedef CGAL::Mesh_criteria_3<Del> Tets_criteria;
typedef CGAL::Implicit_surfaces_mesher_3<Del, Oracle,
                                         Criteria,
                                         Tets_criteria> Mesher;

/////////////// Main function /////////////// 

int main(int argc, char **argv) {

  // Function
  Func F;

  // Oracle (NB: parity oracle is toggled)
  Oracle O (F, K::Point_3 (0,0,0), 4, 1e-3, true);  

  // 2D-complex in 3D-Delaunay triangulation
  Del T;

  // Initial point sample
  Oracle::Points initial_point_sample = O.random_points (10);
  T.insert (initial_point_sample.begin(), initial_point_sample.end());
  
  // Meshing criteria
  CGAL::Chew_4_surfaces::Curvature_size_criterion<Del> c_s_crit (0.1);
  CGAL::Chew_4_surfaces::Uniform_size_criterion<Del> u_s_crit (0.3);
  CGAL::Chew_4_surfaces::Aspect_ratio_criterion<Del> a_r_crit (30);
  std::vector<Criterion*> crit_vect(2);
  crit_vect[0] = &u_s_crit;
  crit_vect[1] = &a_r_crit;
  Criteria C (crit_vect);

  Tets_criteria tets_criteria(1.0);

  std::cout << "Initial number of points: " << T.number_of_vertices() 
            << std::endl;
  
  // Surface meshing
  Mesher mesher (T, O, C, tets_criteria);
//   Chew mesher (T, O, a_r_crit);
  mesher.refine_mesh();
  
  std::cout << "Final number of points: " << T.number_of_vertices() 
            << std::endl;

  // Output
  if (argc >= 2) {
    std::cout << "Writing output to file " << argv[1] << "...";
    std::ofstream os(argv[1]);
    output_surface_facets_to_off (os, T);
    os.close();
    std::cout << " done\n";
  }
}
