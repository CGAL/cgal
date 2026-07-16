#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/AABB_trees/intersection.h>
#include <CGAL/IO/polygon_mesh_io.h>

using Epick = CGAL::Exact_predicates_inexact_constructions_kernel;
using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;
using SCD = CGAL::Simple_cartesian<double>;

template<typename K>
void test()
{
  using P = typename K::Point_3;
  using V = typename K::Vector_3;
  using M =  CGAL::Surface_mesh<P>;
  using fd = typename boost::graph_traits<M>::face_descriptor;
  using Primitive = CGAL::AABB_face_graph_triangle_primitive<M>;
  using Traits = CGAL::AABB_traits_3<K, Primitive>;
  using Tree = CGAL::AABB_tree<Traits>;
  using Aff_tr = CGAL::Aff_transformation_3<K>;

  M m1, m2;
  if(!CGAL::IO::read_polygon_mesh(CGAL::data_file_path("meshes/knot1.off"), m1)){
    std::cout << "error reading knot1" << std::endl;
    exit(1);
  }
  if(!CGAL::IO::read_polygon_mesh(CGAL::data_file_path("meshes/lion.off"), m2)){
    std::cout << "error reading lion" << std::endl;
    exit(1);
  }

  Tree tree1(faces(m1).first, faces(m1).second, m1);
  Tree tree2(faces(m2).first, faces(m2).second, m2);
  tree1.build();
  tree2.build();

  std::vector< std::pair<fd, fd> > inter;
  assert(CGAL::AABB_trees::do_intersect(tree1, tree2));

  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(tree1, tree2, std::back_inserter(inter));
  assert(inter.size() == 1191);
  inter.clear();

  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(tree1, tree2, std::back_inserter(inter), CGAL::parameters::transformation(Aff_tr(CGAL::Translation(), V(1,0,0))));
  assert(inter.size() == 0);
  inter.clear();

  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(tree1, tree2, std::back_inserter(inter), CGAL::parameters::transformation(Aff_tr(0, 1, 0, 1, 0, 0, 0, 0, 1, 1)));
  assert(inter.size() == 1289);
}

#ifdef CGAL_LINKED_WITH_TBB
template<typename K>
void test_parallel()
{
  using P = typename K::Point_3;
  using V = typename K::Vector_3;
  using M =  CGAL::Surface_mesh<P>;
  using fd = typename boost::graph_traits<M>::face_descriptor;
  using Primitive = CGAL::AABB_face_graph_triangle_primitive<M>;
  using Traits = CGAL::AABB_traits_3<K, Primitive>;
  using Tree = CGAL::AABB_tree<Traits>;
  using Aff_tr = CGAL::Aff_transformation_3<K>;

  M m1, m2;
  if(!CGAL::IO::read_polygon_mesh(CGAL::data_file_path("meshes/knot1.off"), m1)){
    std::cout << "error reading knot1" << std::endl;
    exit(1);
  }
  if(!CGAL::IO::read_polygon_mesh(CGAL::data_file_path("meshes/lion.off"), m2)){
    std::cout << "error reading lion" << std::endl;
    exit(1);
  }

  Tree tree1(faces(m1).first, faces(m1).second, m1);
  Tree tree2(faces(m2).first, faces(m2).second, m2);
  tree1.template build<CGAL::Parallel_tag>();
  tree2.template build<CGAL::Parallel_tag>();

  tbb::concurrent_vector< std::pair<fd, fd> > inter;
  assert(CGAL::AABB_trees::do_intersect(tree1, tree2, CGAL::parameters::concurrency_tag(CGAL::Parallel_tag())));

  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(tree1, tree2, std::back_inserter(inter), CGAL::parameters::concurrency_tag(CGAL::Parallel_tag()));
  std::cout << inter.size() << std::endl;
  assert(inter.size() == 1191);
  inter.clear();

  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(tree1, tree2, std::back_inserter(inter), CGAL::parameters::concurrency_tag(CGAL::Parallel_tag()).transformation(Aff_tr(CGAL::Translation(), V(0.5,0,0))));
  std::cout << inter.size() << std::endl;
  assert(inter.size() == 0);
  inter.clear();

  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(tree1, tree2, std::back_inserter(inter), CGAL::parameters::concurrency_tag(CGAL::Parallel_tag()).transformation(Aff_tr(0, 1, 0, 1, 0, 0, 0, 0, 1, 1)));
  std::cout << inter.size() << std::endl;
  assert(inter.size() == 1289);
}
#endif

int main(){
  test<Epick>();
  test<Epeck>();
  test<SCD>();

  return 0;
}
