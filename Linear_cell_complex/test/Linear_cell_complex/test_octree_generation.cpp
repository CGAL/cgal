#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_octree_generation.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <cassert>

void trace_test_begin()
{
  static unsigned int nbtest = 0;
  std::cout<<"Test "<<nbtest++<<" ..."<<std::flush;
}

void trace_test_end()
{
  std::cout<<"Ok."<<std::endl;
}

void trace_display_msg(const char* msg)
{
  std::cout<<"***************** "<<msg<<"***************** "<<std::endl;
}

template<typename LCC>
bool check_octree_basic_properties(LCC& lcc, unsigned int min_darts)
{
  if (!lcc.is_valid())
  {
    std::cout<<"ERROR: the lcc is not valid."<<std::endl;
    return false;
  }

  if (lcc.number_of_darts() < min_darts)
  {
    std::cout<<"ERROR: expected at least "<<min_darts<<" darts, got "
             <<lcc.number_of_darts()<<std::endl;
    return false;
  }

 // Verify that all volumes are hexahedra (6 faces)
  for (auto it = lcc.template one_dart_per_cell<3>().begin();
       it != lcc.template one_dart_per_cell<3>().end(); ++it)
  {
    unsigned int face_count = 0;
    for (auto it2 = lcc.template one_dart_per_incident_cell<2,3>(it).begin();
         it2 != lcc.template one_dart_per_incident_cell<2,3>(it).end(); ++it2)
    {
      ++face_count;
    }

    if (face_count != 6)
    {
      std::cout<<"ERROR: found volume with "<<face_count<<" faces (expected 6)"<<std::endl;
      return false;
    }
  }

  return true;
}

template<typename LCC>
bool test_octree_generation()
{
  trace_test_begin();
  LCC lcc;


  CGAL::compute_octree(lcc, "data/meshes/cube.off");

  if (!check_octree_basic_properties(lcc, 48)) // 2x2x2 = 8 hexaedres × 6 faces × 4 darts
    return false;

  trace_test_end();

  trace_test_begin();
  LCC lcc2;

  // Test with different parameters
  CGAL::compute_octree(lcc2, "data/meshes/cube.off", 3, 3, true);

  if (!check_octree_basic_properties(lcc2, 48))
    return false;

  // Finer grid should produce more voxels
  if (lcc2.number_of_darts() < lcc.number_of_darts())
  {
    std::cout<<"ERROR: finer grid should produce more darts"<<std::endl;
    return false;
  }

  trace_test_end();

  trace_test_begin();
  LCC lcc3;

  // Test with create_all_voxels
  CGAL::compute_octree(lcc3, "data/meshes/tetrahedron.off", 2, 3, true);

  if (!check_octree_basic_properties(lcc3, 24))
    return false;

  trace_test_end();

  trace_test_begin();
  LCC lcc4;

  // Test error handling
  CGAL::compute_octree(lcc4, "non_existent_file.off");

  if (lcc4.number_of_darts() != 0)
  {
    std::cout<<"ERROR: should have 0 darts with invalid file"<<std::endl;
    return false;
  }

  // ****************** TEST SUBDIVISION LEVELS ******************
  trace_test_begin();
  LCC lcc_level1, lcc_level2;

  // Signature: (lcc, filename, initial_grid_size, max_subdivision_level, create_all_voxels, no_remove_outside, regularized)
  // Use different initial_grid_size values to force different dart counts.
  CGAL::compute_octree(lcc_level1, "data/meshes/cube.off", 2, 1, false); // initial_grid_size = 2
  CGAL::compute_octree(lcc_level2, "data/meshes/cube.off", 3, 1, false); // initial_grid_size = 3

  if (!check_octree_basic_properties(lcc_level1, 24)) return false;
  if (!check_octree_basic_properties(lcc_level2, 24)) return false;

  if (lcc_level2.number_of_darts() <= lcc_level1.number_of_darts()) {
    std::cout<<"ERROR: larger initial_grid_size should yield >= darts"<<std::endl;
    return false;
  }

  std::cout<<"Grid size test: "<<lcc_level1.number_of_darts()
           <<" vs "<<lcc_level2.number_of_darts()<<" darts"<<std::endl;
  trace_test_end();

  // ****************** TEST REGULARIZED FLAG  ******************
  trace_test_begin();
  LCC lcc_reg_false, lcc_reg_true;
  CGAL::compute_octree(lcc_reg_false, "data/meshes/cube.off", 2, 1, false, false, false);
  CGAL::compute_octree(lcc_reg_true,  "data/meshes/cube.off", 2, 1, false, false, true);
  if (lcc_reg_false.number_of_darts() != lcc_reg_true.number_of_darts()) {
    std::cout<<"ERROR: regularized stub should not change dart count"<<std::endl;
    return false;
  }
  trace_test_end();

  // ****************** TEST CREATE_ALL_VOXELS EFFECT ******************
  trace_test_begin();
  LCC lcc_inside, lcc_all;
  CGAL::compute_octree(lcc_inside, "data/meshes/sphere.off", 3, 2, false, false, false); // only intersecting
  CGAL::compute_octree(lcc_all,    "data/meshes/sphere.off", 3, 2, true,  false, false); // all voxels

  if (!check_octree_basic_properties(lcc_inside, 24)) return false;
  if (!check_octree_basic_properties(lcc_all, 24))   return false;
  if (lcc_all.number_of_darts() < lcc_inside.number_of_darts()) {
    std::cout<<"ERROR: create_all_voxels should generate more or equal darts ("
             <<lcc_inside.number_of_darts()<<" vs "<<lcc_all.number_of_darts()<<")"<<std::endl;
    return false;
  }
  trace_test_end();

  // ****************** TEST FACEGRAPH OVERLOAD ******************
  trace_test_begin();
  {
    std::ifstream in_fg("data/meshes/cube.off");
    if(!in_fg) { std::cout<<"ERROR: cannot open cube.off\n"; return false; }
    CGAL::Surface_mesh< CGAL::Simple_cartesian<double>::Point_3 > sm;
    if(!(in_fg >> sm)) { std::cout<<"ERROR: invalid cube.off\n"; return false; }
    if(!CGAL::is_triangle_mesh(sm))
      CGAL::Polygon_mesh_processing::triangulate_faces(sm);

    LCC lcc_fg;
    CGAL::compute_octree(lcc_fg, sm, 2, 1, false);
    if(!check_octree_basic_properties(lcc_fg, 48)) return false;
  }
  trace_test_end();

  return true;
}

int main()
{
  std::cout<<"test_octree_generation."<<std::flush;

  // ****************** TEST FOR CMAP ******************
  trace_display_msg("test_octree_generation<LCC3>");
  typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC3;
  if (!test_octree_generation<LCC3>())
  {
    std::cout<<" Error during test_octree_generation<LCC3>."<<std::endl;
    return EXIT_FAILURE;
  }

  // ****************** TEST FOR GMAP ******************
  trace_display_msg("test_octree_generation<GLCC3>");
  typedef CGAL::Linear_cell_complex_for_generalized_map<3> GLCC3;
  if (!test_octree_generation<GLCC3>())
  {
    std::cout<<" Error during test_octree_generation<GLCC3>."<<std::endl;
    return EXIT_FAILURE;
  }

  std::cout<<" Success."<<std::endl;
  return EXIT_SUCCESS;
}