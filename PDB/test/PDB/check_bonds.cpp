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
  CGAL_PDB_FOREACH(CGAL::PDB::Chain::Atom_key k,
                   CGAL::PDB::make_key_range(helix.atoms())) {
    std::cout << k.atom_key() << " " << k.monomer_key() << std::endl;
  }
  CGAL_PDB_FOREACH(CGAL::PDB::Chain::Bond b,
                   helix.bonds()) {
    std::cout << b.first.key().atom_key() << " "
              << b.first.key().monomer_key() << ": "
              << b.second.key().atom_key() << " "
              << b.second.key().monomer_key()
              << std::endl;
  }
  /*CGAL_PDB_FOREACH(CGAL::PDB::Atom k,
                   CGAL::PDB::make_atom_range(helix.atoms())) {
    std::cout << ".sphere " << k.point() << " ./15" << std::endl;
  }
  CGAL_PDB_FOREACH(CGAL::PDB::Chain::Bond b,
                   helix.bonds()) {
    std::cout << ".cylinder " << b.first.atom().point() << " "
              << b.second.atom().point() << " .1"
              << std::endl;
              } */
  std::cout << bsz << " " << ptsz << std::endl;
  assert(bsz == ptsz-1);
  std::cout << bsz << " " << ptsz << std::endl;
  return EXIT_SUCCESS;
}
