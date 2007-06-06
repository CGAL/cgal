#include <CGAL/PDB/Protein.h>
#include <CGAL/PDB/geometry.h>
#include <fstream>
#include <vector>
#include <cassert>

int main(int, char *[]){
  std::ifstream input("data/check_bonds.pdb");
  CGAL_PDB_NS::Protein helix(input);
  //helix.dump(std::cout);
  std::vector<CGAL_PDB_NS::Point> points;
  std::vector<CGAL_PDB_NS::Protein::Bond> bonds;
  CGAL_PDB_NS::backbone_coordinates_and_bonds(helix, std::back_inserter(points), std::back_inserter(bonds));
  assert(bonds.size() == points.size()-2);
  return EXIT_SUCCESS;
}
