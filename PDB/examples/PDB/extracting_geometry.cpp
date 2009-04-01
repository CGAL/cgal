#include <CGAL/PDB/PDB.h>
#include <CGAL/PDB/range.h>
#include <fstream>

template <class Vector, class Range>
void insert(Vector &v, Range r) {
  v.insert(v.end(), r.begin(), r.end());
}

int main(int, char *[]) {
  std::ifstream in("data/simple.pdb");
  CGAL::PDB::PDB pdb(in);
  CGAL::PDB::Model m=pdb.models().begin()->model();
  CGAL::PDB::Chain c=m.chains().begin()->chain();
  // all of these can be run on the Model or on the Chain (or Monomer)
  {
    // get all spheres
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    K k;
    std::vector<K::Sphere_3> points;
    insert(points, CGAL::PDB::make_sphere_3_range(CGAL::PDB::make_atom_range(m.atoms()), k));
    std::cout << points.size() << std::endl;
  }
  {
    // get backbone points
    std::vector<CGAL::PDB::Point> points;
    insert(points,
           CGAL::PDB::make_point_range(CGAL::PDB::make_atom_range(CGAL::PDB::make_backbone_range(m.atoms())))
           );
    std::cout << points.size() << std::endl;
  }
  {
    // get all points and bonds
    CGAL::PDB::index_atoms(c);
    std::vector<CGAL::PDB::Point> points;
    std::vector<std::pair<unsigned int, unsigned int> > bonds;
    insert(points,
           CGAL::PDB::make_point_range(CGAL::PDB::make_atom_range(c.atoms())));
    insert(bonds,
           CGAL::PDB::make_bond_indices_range(c.bonds()));
    for (unsigned int i=0; i< bonds.size(); ++i) {
      CGAL_assertion(bonds[i].first < points.size());
      CGAL_assertion(bonds[i].second < points.size());
    }
    std::cout << points.size() << " " << bonds.size() << std::endl;
  }

  return EXIT_SUCCESS;
}
