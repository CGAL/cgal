#include <CGAL/Rigid_triangle_mesh_collision_detection.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Surface_mesh;
typedef CGAL::Polyhedron_3<K>                      Polyhedron_3;
namespace params = CGAL::parameters;

void test_remove()
{
  std::cout << "test_remove()" << std::endl;
  Surface_mesh sm;
  CGAL::Rigid_triangle_mesh_collision_detection<Surface_mesh> collision_detection;
  collision_detection.add_mesh(sm); // 0
  collision_detection.add_mesh(sm); // 1
  collision_detection.add_mesh(sm); // 2
  collision_detection.add_mesh(sm); // 3
  collision_detection.add_mesh(sm); // 4

  assert( collision_detection.is_valid_index(0) );
  assert( collision_detection.is_valid_index(1) );
  assert( collision_detection.is_valid_index(2) );
  assert( collision_detection.is_valid_index(3) );
  assert( collision_detection.is_valid_index(4) );
  assert(!collision_detection.is_valid_index(5) );

  assert(collision_detection.size() == 5);

  collision_detection.remove_mesh(4);
  assert(collision_detection.size() == 4);
  assert(!collision_detection.is_valid_index(4));

  std::size_t id = collision_detection.add_mesh(sm);
  assert(id == 4);

  collision_detection.remove_mesh(2);
  assert(collision_detection.size() == 4);
  assert(!collision_detection.is_valid_index(2));
  id = collision_detection.add_mesh(sm);
  assert(id == 2);

  assert( collision_detection.is_valid_index(0) );
  assert( collision_detection.is_valid_index(1) );
  assert( collision_detection.is_valid_index(2) );
  assert( collision_detection.is_valid_index(3) );
  assert( collision_detection.is_valid_index(4) );

  collision_detection.remove_mesh(3);
  assert( collision_detection.size()==4);
  assert( collision_detection.size_of_garbage()==1);
  collision_detection.remove_mesh(3);
  assert( collision_detection.size()==4);
  assert( collision_detection.size_of_garbage()==1);
  collision_detection.remove_mesh(2);
  assert( collision_detection.size()==3);
  assert( collision_detection.size_of_garbage()==2);
  collision_detection.remove_mesh(1);
  assert( collision_detection.size()==2);
  assert( collision_detection.size_of_garbage()==3);
  collision_detection.remove_mesh(4);
  assert( collision_detection.size()==1);
  assert( collision_detection.size_of_garbage()==4);
  collision_detection.remove_mesh(0);
  assert( collision_detection.size()==0);
  assert( collision_detection.size_of_garbage()==5);
}

template <class TriangleMesh, class Index>
void test_intersections(Index index, const char* type)
{
  std::cout << "test_intersections<"<<type<<">()" << std::endl;
  TriangleMesh tm1, tm2, tm3;
  std::ifstream input("data/small_spheres.off");
  assert(input);
  input >> tm1;
  input.close();
  input.open(CGAL::data_file_path("meshes/blobby.off"));
  assert(input);
  input >> tm2;
  input.close();
  input.open("data-coref/large_cube_coplanar.off");
  assert(input);
  input >> tm3;
  input.close();

  CGAL::Rigid_triangle_mesh_collision_detection<TriangleMesh> collision_detection;
  collision_detection.add_mesh(tm1, params::face_index_map(get(index, tm1))); // 0
  // add tm1 using an external tree
  typename CGAL::Rigid_triangle_mesh_collision_detection<TriangleMesh>::AABB_tree
    tm1_tree(std::begin(faces(tm1)), std::end(faces(tm1)), tm1);
  collision_detection.add_mesh(tm1_tree, tm1, params::face_index_map(get(index, tm1))); // 1 small_spheres
  collision_detection.add_mesh(tm2, params::face_index_map(get(index, tm2))); // 2 blobby
  collision_detection.add_mesh(tm3, params::face_index_map(get(index, tm3))); // 3 large_cube_coplanar
  // pool is 0 1 2 3
  collision_detection.remove_mesh(0);
  // pool is now 3 1 2
  std::vector<std::size_t> inter_res;
  std::vector< std::pair<std::size_t, bool> > inter_and_inclus_res;

  // spheres intersects both cube and blobby
  inter_res = collision_detection.get_all_intersections(1);
  assert(inter_res.size() == 2);
  assert(inter_res[0] == 3);
  assert(inter_res[1] == 2);

  // blobby intersects only spheres
  inter_res = collision_detection.get_all_intersections(2);
  assert(inter_res.size() == 1);
  assert(inter_res[0] == 1);

  // blobby is included into cube
  inter_and_inclus_res = collision_detection.get_all_intersections_and_inclusions(2);
  assert(inter_and_inclus_res.size() == 2);
  assert(inter_and_inclus_res[0].first == 3 && inter_and_inclus_res[0].second);
  assert(inter_and_inclus_res[1].first == 1 && !inter_and_inclus_res[1].second);

  // set transformations
  K::Aff_transformation_3 tm1_transf(-0.22115682390903474, 0.97263618730220114, 0.071193443439028392,
                                     -2.3805504333121519, 0.37717228816945586, 0.017983434697950496,
                                     0.92596849898551192, -4.6543355030238258, 0.89935016777420218,
                                     0.23163644624001295, -0.37082857562196114, -5.4867495857613582);

  K::Aff_transformation_3 tm2_transf(-0.10062050938738221, -0.92610032901754502, -0.36361200981845854,
                                     -2.1218726603809985, 0.98718082597554757, -0.13844084127856282,
                                     0.079423864753106477, -4.6627185424110094, -0.1238932198179533,
                                     -0.35095913445824267, 0.92815847570523291, -5.349009707055858);

  K::Aff_transformation_3 tm3_transf(0.8523140340417823, -0.41651563577194606, -0.3163471392835962,
                                     -2.0549216027623891, 0.34883891716890358, 0.90335137983769587,
                                     -0.24953495629622252, -5.8394794922248607, 0.3897078357485903,
                                     0.10232795171810326, 0.91523592207328008, -5.2436990668816268);

  collision_detection.set_transformation(1, tm1_transf);
  collision_detection.set_transformation(2, tm2_transf);
  collision_detection.set_transformation(3, tm3_transf);

  // spheres intersects blobby
  inter_res = collision_detection.get_all_intersections(1);
  assert(inter_res.size() == 1);
  assert(inter_res[0] == 2);

  // blobby intersect spheres and cube
  inter_res = collision_detection.get_all_intersections(2);
  assert(inter_res.size() == 2);
  assert(inter_res[0] == 3);
  assert(inter_res[1] == 1);

  // cube contains one CC of spheres
  inter_and_inclus_res = collision_detection.get_all_intersections_and_inclusions(3);
  assert(inter_and_inclus_res.size() == 2);
  assert(inter_and_inclus_res[0].first == 1 && inter_and_inclus_res[0].second);
  assert(inter_and_inclus_res[1].first == 2 && !inter_and_inclus_res[1].second);

  // one CC of spheres is contained by cube
  inter_and_inclus_res = collision_detection.get_all_intersections_and_inclusions(1);
  assert(inter_and_inclus_res.size() == 2);
  assert(inter_and_inclus_res[0].first == 3 && inter_and_inclus_res[0].second);
  assert(inter_and_inclus_res[1].first == 2 && !inter_and_inclus_res[1].second);
}

int main()
{
  test_remove();
  test_intersections<Surface_mesh>(boost::face_index, "Surface_mesh");
  test_intersections<Polyhedron_3>(boost::face_external_index, "Polyhedron_3");

  return 0;
}
