
#include <CGAL/Random.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#include <iostream>
#include <fstream>


template<typename T3>
bool load_binary_triangulation(std::istream& is, T3& t3)
{
  std::string s;
  if (!(is >> s)) return false;
  bool binary = (s == "binary");
  if (binary) {
    if (!(is >> s)) return false;
  }
  if (s != "CGAL" || !(is >> s) || s != "c3t3")
    return false;

  std::getline(is, s);
  if (binary) CGAL::set_binary_mode(is);
  is >> t3;
  return bool(is);
}

template<typename T3>
bool save_binary_triangulation(std::ostream& os, const T3& t3)
{
  os << "binary CGAL c3t3\n";
  CGAL::set_binary_mode(os);
  return !!(os << t3);
}

template<typename T3>
void save_ascii_triangulation(const char* filename, const T3& t3)
{
  if (!t3.is_valid(true))
    std::cerr << "Invalid triangulation!" << std::endl;

  CGAL::Tetrahedral_remeshing::debug::dump_triangulation_cells(
        t3, filename);
}

