#include <CGAL/PDB/transforms.h>
#include <CGAL/PDB/PDB.h>
#include <CGAL/PDB/Protein.h>
#include <CGAL/PDB/Transform.h>
CGAL_PDB_BEGIN_NAMESPACE

  void transform_protein(const Transform &t, Protein &p) {
    for (Protein::Residues_iterator rit= p.residues_begin(); rit != p.residues_end(); ++rit){
      for (Residue::Atoms_iterator ait= rit->atoms_begin(); ait != rit->atoms_end(); ++ait){
	ait->second.set_cartesian_coords(t(ait->second.cartesian_coords()));
      }
    }
  }


  void transform_pdb(const Transform &t, PDB &pdb) {
    for (unsigned int i=0;i< pdb.number_of_models(); ++i){
      Model &m= pdb.model(i);
      std::cout << "Model " << i << " has " << m.number_of_chains() << " chains."<< std::endl;
      for (unsigned int j=0; j< m.number_of_chains(); ++j){
	Protein &p= m.chain(j);
	transform_protein(t, p);
      }
    }
  }
CGAL_PDB_END_NAMESPACE
