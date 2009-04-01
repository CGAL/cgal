#include <CGAL/PDB/Chain.h>
#include <CGAL/PDB/PDB.h>
#include <CGAL/PDB/range.h>
#include <fstream>
#include <vector>
#include <cassert>

int main(int, char *[]){
  std::ifstream input("data/check_bonds.pdb");
  CGAL::PDB::PDB pdb(input);
  CGAL::PDB::Chain helix= pdb.models().begin()->model().chains().begin()->chain();
  //helix.dump(std::cout);
  std::vector<CGAL::PDB::Point> points;
  std::vector<CGAL::PDB::Chain::Bond> bonds;
  unsigned int ptsz=CGAL::PDB::size(helix.atoms());
  unsigned int bsz=CGAL::PDB::size(helix.bonds());
  assert(bsz == ptsz-2);
  std::cout << bsz << " " << ptsz << std::endl;
  return EXIT_SUCCESS;
}
