
#include <iostream>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_traits_construct_by_sorting.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Timer.h>

static std::size_t R = 100;

template<typename K>
double benchmark_recursive_partitioning_construction(std::string input_path) {
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  std::ifstream in(input_path);
  Polyhedron polyhedron;
  in >> polyhedron;

  CGAL::Timer t;
  t.start();
  for (int i = 0; i < R; ++i) {
    Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
    tree.build();
  }
  t.stop();

  return t.time() / R;
}

template<typename K>
double benchmark_recursive_partitioning_traversal(std::string input_path) {

  return 0;
}

template<typename K>
double benchmark_sorting_construction(std::string input_path) {
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
  typedef CGAL::AABB_traits_construct_by_sorting<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  std::ifstream in(input_path);
  Polyhedron polyhedron;
  in >> polyhedron;

  CGAL::Timer t;
  t.start();
  for (int i = 0; i < R; ++i) {
    Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
    tree.build();
  }
  t.stop();

  return t.time() / R;
}

template<typename K>
double benchmark_sorting_traversal(std::string input_path) {

  return 0;
}

int main(int argc, char **argv) {

  // Determine our data source, with a default if no path is provided
  std::string input_path = argc > 1 ? argv[1] : "data/handle.off";

  std::cout << "{| class=\"wikitable\"\n";
  std::cout << "! Simple Cartesian float !! Recursive Partition !! Hilbert Sort" << "\n|-\n";
  std::cout << "| construction || " << std::flush
            << benchmark_recursive_partitioning_construction<CGAL::Simple_cartesian<float>>(input_path)
            << " s || " << std::flush
            << benchmark_sorting_construction<CGAL::Simple_cartesian<float>>(input_path)
            << " s" << "\n|-\n";
  std::cout << "| traversal || " << std::flush
            << benchmark_recursive_partitioning_traversal<CGAL::Simple_cartesian<float>>(input_path)
            << " s || " << std::flush
            << benchmark_sorting_traversal<CGAL::Simple_cartesian<float>>(input_path)
            << " s" << "\n|-\n";
  std::cout << "|}\n\n";


  std::cout << "{| class=\"wikitable\"\n";
  std::cout << "! Cartesian float !! Recursive Partition !! Hilbert Sort" << "\n|-\n";
  std::cout << "| construction || " << std::flush
            << benchmark_recursive_partitioning_construction<CGAL::Cartesian<float>>(input_path)
            << " s || " << std::flush
            << benchmark_sorting_construction<CGAL::Cartesian<float>>(input_path)
            << " s" << "\n|-\n";
  std::cout << "| traversal || " << std::flush
            << benchmark_recursive_partitioning_traversal<CGAL::Cartesian<float>>(input_path)
            << " s || " << std::flush
            << benchmark_sorting_traversal<CGAL::Cartesian<float>>(input_path)
            << " s" << "\n|-\n";
  std::cout << "|}\n\n";


  std::cout << "{| class=\"wikitable\"\n";
  std::cout << "! Simple Cartesian double !! Recursive Partition !! Hilbert Sort" << "\n|-\n";
  std::cout << "| construction || " << std::flush
            << benchmark_recursive_partitioning_construction<CGAL::Simple_cartesian<double>>(input_path)
            << " s || " << std::flush
            << benchmark_sorting_construction<CGAL::Simple_cartesian<double>>(input_path)
            << " s" << "\n|-\n";
  std::cout << "| traversal || " << std::flush
            << benchmark_recursive_partitioning_traversal<CGAL::Simple_cartesian<double>>(input_path)
            << " s || " << std::flush
            << benchmark_sorting_traversal<CGAL::Simple_cartesian<double>>(input_path)
            << " s" << "\n|-\n";
  std::cout << "|}\n\n";


  std::cout << "{| class=\"wikitable\"\n";
  std::cout << "! Cartesian double !! Recursive Partition !! Hilbert Sort" << "\n|-\n";
  std::cout << "| construction || " << std::flush
            << benchmark_recursive_partitioning_construction<CGAL::Cartesian<double>>(input_path)
            << " s || " << std::flush
            << benchmark_sorting_construction<CGAL::Cartesian<double>>(input_path)
            << " s" << "\n|-\n";
  std::cout << "| traversal || " << std::flush
            << benchmark_recursive_partitioning_traversal<CGAL::Cartesian<double>>(input_path)
            << " s || " << std::flush
            << benchmark_sorting_traversal<CGAL::Cartesian<double>>(input_path)
            << " s" << "\n|-\n";
  std::cout << "|}\n\n";


  std::cout << "{| class=\"wikitable\"\n";
  std::cout << "! Epic !! Recursive Partition !! Hilbert Sort" << "\n|-\n";
  std::cout << "| construction || " << std::flush
            << benchmark_recursive_partitioning_construction<CGAL::Exact_predicates_inexact_constructions_kernel>(
                    input_path)
            << " s || " << std::flush
            << benchmark_sorting_construction<CGAL::Exact_predicates_inexact_constructions_kernel>(input_path)
            << " s" << "\n|-\n";
  std::cout << "| traversal || " << std::flush

            << benchmark_recursive_partitioning_traversal<CGAL::Exact_predicates_inexact_constructions_kernel>(
                    input_path)
            << " s || " << std::flush
            << benchmark_sorting_traversal<CGAL::Exact_predicates_inexact_constructions_kernel>(input_path)
            << " s" << "\n|-\n";
  std::cout << "|}\n\n";


}
