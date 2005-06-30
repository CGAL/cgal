#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Complex_2_in_triangulation_surface_mesh_cell_base_3.h>
#include <CGAL/Mesh_3/Complex_2_in_triangulation_cell_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Implicit_surfaces_mesher_3.h>
#include <CGAL/Surface_mesher/Criteria/Standard_criteria.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <CGAL/IO/File_medit.h>

#include <CGAL/Surface_mesher/Oracles/Implicit_oracle.h>

#include <CGAL/Mesh_criteria_3.h>

#include "debug.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "Isosurface.h"
#include "parameters.h"

/////////////// Types /////////////// 

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Triangulation_vertex_base_with_info_3<bool, K> Vb;
typedef CGAL::Complex_2_in_triangulation_surface_mesh_cell_base_3<K> CCb;
typedef CGAL::Mesh_3::Complex_2_in_triangulation_cell_base_3<K, CCb> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Del;
typedef CGAL::Inrimage_isosurface<K::FT> Isosurface;
typedef CGAL::Surface_mesher::Implicit_oracle<K, Isosurface> Oracle;
typedef CGAL::Surface_mesher::Refine_criterion<Del> Criterion;
typedef CGAL::Surface_mesher::Standard_criteria <Criterion > Criteria;
typedef CGAL::Mesh_criteria_3<Del> Tets_criteria;
typedef CGAL::Implicit_surfaces_mesher_3<Del, Oracle,
                                         Criteria,
                                         Tets_criteria> Mesher;

/////////////// Main function /////////////// 

void usage(char * name, std::string error = std::string() )
{
  if( ! error.empty() )
    std::cout << "Error: " << error << "\n";

  std::cout << 
    "Usage:\n"
    "  " << name << " inrimage_file iso_value [out_file.mesh]\n";
  exit(1);
}

int main(int argc, char **argv) {
  if(argc <=1 )
    usage(argv[0]);
  
  std::stringstream argv2(argv[2]);
  int iso_value;
  if( ! (argv2 >> iso_value) )
    usage(argv[0], "bad iso_value");
  Isosurface surface(argv[1], iso_value);

  // Oracle (NB: parity oracle is toggled)
  Oracle O (surface,
	    K::Point_3 (0,0,0),
	    enclosing_sphere_radius,
	    precision,
	    bipolar_oracle);

  // 2D-complex in 3D-Delaunay triangulation
  Del T;

  // Initial point sample
  Oracle::Points initial_point_sample = 
    O.random_points (number_of_initial_points);
  T.insert (initial_point_sample.begin(), initial_point_sample.end());
  
  // Meshing criteria
  CGAL::Surface_mesher::Curvature_size_criterion<Del> 
    c_s_crit (curvature_bound);
//   CGAL::Surface_mesher::Uniform_size_criterion<Del>
//     u_s_crit (size_bound);
  CGAL::Surface_mesher::Aspect_ratio_criterion<Del>
    a_r_crit (aspect_ratio_bound);

  std::vector<Criterion*> crit_vect;
  crit_vect.push_back (&c_s_crit);
//   crit_vect.push_back (&u_s_crit);
//   crit_vect.push_back(&a_r_crit);
  Criteria C (crit_vect);

  Tets_criteria tets_criteria(tets_aspect_ratio_bound, tets_size_bound);

  std::cout << "Initial number of points: " << T.number_of_vertices() 
            << std::endl;
  
  // Surface meshing
  Mesher mesher (T, O, C, tets_criteria);
  mesher.refine_surface();
  mesher.refine_mesh();
//   int i = 100;
//   while(!mesher.done())
//     {
//       std::stringstream s;
//       s << "out." << i++ << ".mesh";
//       std::cerr << s.str() << std::endl;
//       std::ofstream os3(s.str().c_str());
//       output_pslg_to_medit(os3, T);
//       mesher.step_by_step();
//     }
//   exit(0);
  std::cout << "Final number of points: " << T.number_of_vertices() 
            << std::endl;

  // Output
  if (argc >= 3) {
    std::ofstream os(argv[2]);
    output_pslg_to_medit(os, T);
    os.close();
  }
}
