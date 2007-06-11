#include <CGAL/PDB/PDB.h>
#include <CGAL/PDB/iterator.h>
#include <fstream>



int main(int, char *[]) {
  std::ifstream in("data/simple.pdb");
  CGAL::PDB::PDB pdb(in);
  CGAL::PDB::Model m=pdb.models_begin()->model();
  CGAL::PDB::Chain c=m.chains_begin()->chain();
  // all of these can be run on the Model or on the Chain (or Monomer)
  {
    // get all spheres
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    K k;
    std::vector<K::Sphere_3> points;
    points.insert(points.end(),
		  CGAL::PDB::make_sphere_3_iterator(CGAL::PDB::make_atom_iterator(m.atoms_begin()), k),
		  CGAL::PDB::make_sphere_3_iterator(CGAL::PDB::make_atom_iterator(m.atoms_end()), k)
		  );
    std::cout << points.size() << std::endl;
  }
  {
    // get backbone points
    std::vector<CGAL::PDB::Point> points;
    points.insert(points.end(),
		  CGAL::PDB::make_point_iterator(CGAL::PDB::make_atom_iterator(CGAL::PDB::make_backbone_iterator(m.atoms_begin(), m.atoms_end()))),
		  CGAL::PDB::make_point_iterator(CGAL::PDB::make_atom_iterator(CGAL::PDB::make_backbone_iterator(m.atoms_end(), m.atoms_end())))
		  );
    std::cout << points.size() << std::endl;
  }
  {
    // get all points and bonds
    CGAL::PDB::index_atoms(c);
    std::vector<CGAL::PDB::Point> points;
    std::vector<std::pair<unsigned int, unsigned int> > bonds;
    points.insert(points.end(),
		  CGAL::PDB::make_point_iterator(CGAL::PDB::make_atom_iterator(c.atoms_begin())),
		  CGAL::PDB::make_point_iterator(CGAL::PDB::make_atom_iterator(c.atoms_end())));
    bonds.insert(bonds.end(),
		 CGAL::PDB::make_bond_indices_iterator(c.bonds_begin()),
		 CGAL::PDB::make_bond_indices_iterator(c.bonds_end()));
    for (unsigned int i=0; i< bonds.size(); ++i) {
      CGAL_assertion(bonds[i].first < points.size());
      CGAL_assertion(bonds[i].second < points.size());
    }
    std::cout << points.size() << " " << bonds.size() << std::endl;
  }

  return EXIT_SUCCESS;
}
