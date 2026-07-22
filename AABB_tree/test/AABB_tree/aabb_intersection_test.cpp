#include <fstream>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_mesh_io.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_trees/intersection.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>

#include <CGAL/Polyhedron_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

void test(const std::string fname1, const std::string fname2, std::size_t nb_inter)
{
  typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
  typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> Primitive;
  typedef CGAL::AABB_traits_3<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  Surface_mesh sm1, sm2;
  CGAL::IO::read_polygon_mesh(fname1, sm1);
  CGAL::IO::read_polygon_mesh(fname2, sm2);

  Tree tree1(faces(sm1).begin(), faces(sm1).end(), sm1);
  Tree tree2(faces(sm2).begin(), faces(sm2).end(), sm2);

  std::vector< std::pair<Primitive::Id, Primitive::Id> > inter;
  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(tree1, tree2, std::back_inserter(inter));
  assert( CGAL::AABB_trees::do_intersect(tree1, tree2) == (nb_inter!=0) );
  assert( inter.size() == nb_inter );
}

int main()
{
    test(CGAL::data_file_path("meshes/cube.off"), CGAL::data_file_path("meshes/cylinder.off"), 121);
    test(CGAL::data_file_path("meshes/cube.off"), CGAL::data_file_path("meshes/femur.off"), 0);
    test(CGAL::data_file_path("meshes/cube.off"), CGAL::data_file_path("meshes/pinion_small.off"), 0);
    test(CGAL::data_file_path("meshes/cylinder.off"), CGAL::data_file_path("meshes/femur.off"), 0);
    test(CGAL::data_file_path("meshes/cylinder.off"), CGAL::data_file_path("meshes/pinion_small.off"), 0);
    test(CGAL::data_file_path("meshes/femur.off"), CGAL::data_file_path("meshes/pinion_small.off"), 905);
    return EXIT_SUCCESS;
}
