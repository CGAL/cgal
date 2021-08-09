
#include <iostream>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

template<typename K>
double benchmark_recursive_partitioning_construction(std::string input_path) {

  return 0;
}

template<typename K>
double benchmark_recursive_partitioning_traversal(std::string input_path) {

  return 0;
}

template<typename K>
double benchmark_sorting_construction(std::string input_path) {

  return 0;
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
