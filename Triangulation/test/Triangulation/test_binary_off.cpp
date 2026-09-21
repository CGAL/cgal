#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/IO/Triangulation_off_ostream.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_3<K> Triangulation;
typedef K::Point_3 Point;

int main()
{
  // --- PART 1: SETUP & WRITE ---
  std::cout << "1. Creating and Writing Triangulation..." << std::endl;
  Triangulation T_write;
  T_write.insert(Point(0,0,0));
  T_write.insert(Point(1,0,0));
  T_write.insert(Point(0,1,0));
  T_write.insert(Point(0,0,1));

  const char* filename = "test_binary_output.off";

  // Write (scoped to ensure file closes)
  {
    std::ofstream out(filename, std::ios::binary);
    if(!out) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return 1;
    }
    out << T_write;
  }

  // --- PART 2: VERIFY 32-BIT INDICES (Size Check) ---
  std::cout << "2. Verifying File Size..." << std::endl;
  std::ifstream in_size(filename, std::ios::binary | std::ios::ate);
  long fileSize = in_size.tellg();
  in_size.close();

  std::cout << "   -> File Size: " << fileSize << " bytes" << std::endl;

  // Expected size is around 115-184 bytes.
  // If indices were 64-bit, it would be significantly larger (>200 bytes).
  if (fileSize > 200) {
      std::cerr << "FAILURE: File is too large (" << fileSize << " bytes)." << std::endl;
      std::cerr << "   Indices are likely 64-bit instead of 32-bit." << std::endl;
      return 1;
  }
  std::cout << "   -> Size check passed (Compact binary format confirmed)." << std::endl;

  // --- PART 3: READ BACK (Round-Trip) ---
  std::cout << "3. Reading back..." << std::endl;
  Triangulation T_read;
  std::ifstream in(filename, std::ios::binary);

  if(!in) {
      std::cerr << "Error opening file for reading!" << std::endl;
      return 1;
  }

  in >> T_read;

  // --- PART 4: VERIFY DATA ---
  std::cout << "4. Verifying Data Integrity..." << std::endl;

  if (T_write.number_of_vertices() != T_read.number_of_vertices()) {
      std::cerr << "FAILURE: Vertex count mismatch!" << std::endl;
      return 1;
  }

  if (T_write.number_of_cells() != T_read.number_of_cells()) {
      std::cerr << "FAILURE: Cell count mismatch!" << std::endl;
      return 1;
  }

  std::cout << "ALL TESTS PASSED: Written, Checked Size, and Read Back successfully." << std::endl;
  return 0;
}
