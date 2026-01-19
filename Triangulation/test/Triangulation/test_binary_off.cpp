#include <iostream>
#include <fstream>
#include <cassert>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/IO/Triangulation_off_ostream.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_3<K> Triangulation;
typedef K::Point_3 Point;

int main()
{
  // --- PART 1: WRITE ---
  std::cout << "1. Creating and Writing Triangulation..." << std::endl;
  Triangulation T_write;
  T_write.insert(Point(0,0,0));
  T_write.insert(Point(1,0,0));
  T_write.insert(Point(0,1,0));
  T_write.insert(Point(0,0,1));

  const char* filename = "test_binary_output.off";
  std::ofstream out(filename, std::ios::binary);

  if(!out) {
      std::cerr << "Error opening file for writing!" << std::endl;
      return 1;
  }

  // Write using the fix
  out << T_write;
  out.close();

  // --- PART 2: READ BACK ---
  std::cout << "2. Reading back..." << std::endl;
  Triangulation T_read;
  std::ifstream in(filename, std::ios::binary);

  if(!in) {
      std::cerr << "Error opening file for reading!" << std::endl;
      return 1;
  }

  // Try to read it back using standard operator
  in >> T_read;

  // --- PART 3: VERIFY ---
  std::cout << "3. Verifying..." << std::endl;
  
  // Check vertex count
  if (T_write.number_of_vertices() != T_read.number_of_vertices()) {
      std::cerr << "FAILURE: Vertex count mismatch!" << std::endl;
      std::cerr << "   Written: " << T_write.number_of_vertices() << std::endl;
      std::cerr << "   Read:    " << T_read.number_of_vertices() << std::endl;
      return 1;
  }

  // Check cell (tetrahedra) count
  if (T_write.number_of_cells() != T_read.number_of_cells()) {
      std::cerr << "FAILURE: Cell count mismatch!" << std::endl;
      return 1;
  }

  std::cout << "VICTORY: Round-trip test passed! (Write -> Read -> Match)" << std::endl;
  return 0;
}
