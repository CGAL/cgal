#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesher/Surface_mesher.h>
// #include <CGAL/Surface_mesher/Surface_mesher_regular_edges.h>
// #include <CGAL/Surface_mesher/Surface_mesher_regular_edges_without_boundary.h>
// #include <CGAL/Surface_mesher/Surface_mesher_manifold.h>
#include <CGAL/Surface_mesher/Criteria/Standard_criteria.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#ifndef SURFACE_MESHER_POLYHEDRAL
#  include <CGAL/Implicit_surface_oracle.h>
#  include "implicit_function.h"
#else // IMPLICIT
#  include <CGAL/Surface_mesher/Oracles/Polyhedral.h>
#endif

#include <fstream>

// triangulation type generator: Kernel -> Tr
templace <class Kernel>
class Surface_mesh_Delaunay_triangulation_generator {
  typedef CGAL::Surface_mesh_vertex_base_3<K> Vb;
  typedef CGAL::Surface_mesh_cell_base_3<K> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
public:
  typedef CGAL::Delaunay_triangulation_3<K, Tds> Type;
  typedef Type type;
};



// Oracle
#ifndef SURFACE_MESHER_POLYHEDRAL
typedef Function <K::FT> Func;
typedef CGAL::Implicit_surface_oracle<K, Func> Oracle;
#else
typedef CGAL::Surface_mesher::Polyhedral <Tr> Oracle;
#endif

typedef CGAL::Surface_mesher::Refine_criterion<Tr> Criterion;
typedef CGAL::Surface_mesher::Standard_criteria <Criterion > Criteria;

typedef CGAL::Surface_mesher::Surface_mesher<Tr, Oracle, Criteria> SM;
// typedef CGAL::Surface_mesher::Surface_mesher_regular_edges<Tr, Oracle, Criteria> SMRE;
// typedef CGAL::Surface_mesher::Surface_mesher_regular_edges_without_boundary<Tr, Oracle, Criteria> SMREWB;
// typedef CGAL::Surface_mesher::Surface_mesher_manifold<Tr, Oracle, Criteria> SMM;
// typedef CGAL::Surface_mesher::Surface_mesher_regular_edges_without_boundary_base<Tr, Oracle, Criteria> SMREWBB;
// typedef CGAL::Surface_mesher::Surface_mesher_manifold<Tr, Oracle, Criteria, 
// CGAL::Surface_mesher::Surface_mesher_manifold_base <Tr, Oracle, Criteria, SMREWBB> 
// > SMMWB;


typedef SM Surface;  // basic mesher
//typedef SMM Surface;  // manifold with boundary
// typedef SMMWB Surface;  // manifold without boundary

// typedef SMRE Surface;  // only regular edges
//typedef SMREWB Surface;  // only regular edges, without boundary


/////////////// Main function ///////////////

int main(int argc, char **argv) {

  Func F;
  Oracle O (F, K::Point_3 (0,0,0), 6, 1e-6, true);  // parity oracle
                                                    // toggled
  // 3D-Delaunay triangulation
  Tr tr;

  // Initial point sample
  O.initial_points (CGAL::inserter(tr), 30);

  // Meshing criteria
  CGAL::Surface_mesher::Curvature_size_criterion<Tr> c_s_crit (10000); 
                                         // bound on Hausdorff distance
                                         // does not play any role if
                                         // bigger than the square of
                                         // the Uniform_size_criterion
  CGAL::Surface_mesher::Uniform_size_criterion<Tr> u_s_crit (0.1); 
                           // bound on radii of surface Delaunay balls
  CGAL::Surface_mesher::Aspect_ratio_criterion<Tr> a_r_crit (30); 
                          // lower bound on minimum angle in degrees
  std::vector<Criterion*> crit_vect;
  crit_vect.push_back (&a_r_crit);
  crit_vect.push_back (&u_s_crit);
  crit_vect.push_back (&c_s_crit);

  Criteria C (crit_vect);

  std::cout << "Initial number of points: " << tr.number_of_vertices() 
            << std::endl;

  // 2D-complex in 3D-Delaunay triangulation
  C2t3 Co2 (tr);
  // Surface meshing
  Surface mesher (tr, Co2, O, C);
//   Surface mesher (tr, O, a_r_crit);
  mesher.refine_mesh(true);
  std::cout << "Final number of points: " << tr.number_of_vertices() 
            << std::endl;

  // Output
#ifdef SURFACE_MESHER_POLYHEDRAL
  const int out_filename_position_in_argv = 2;
#else
  const int out_filename_position_in_argv = 1;
#endif

  if (argc >= out_filename_position_in_argv + 1) {
    std::cout << "Writing output to file " << argv[out_filename_position_in_argv] << "...";
    std::ofstream os(argv[out_filename_position_in_argv]);
    output_surface_facets_to_off (os, tr);
    os.close();
    std::cout << " done\n";
  }

}
