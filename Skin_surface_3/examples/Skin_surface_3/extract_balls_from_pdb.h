#ifndef EXTRACT_BALLS_FROM_PDB_H
#define EXTRACT_BALLS_FROM_PDB_H

#include "dsrpdb/PDB.h"
#include "dsrpdb/geometry.h"

template <class Traits, class OutputIterator>
void extract_balls_from_pdb(char *filename, 
			    Traits const &t,
			    OutputIterator weighted_points) 
{
  std::ifstream in(filename);

  if (!in) {
    std::cerr << "Error opening input file " << filename << std::endl;
  }

  dsrpdb::PDB pdb(in, true);
  dsrpdb::all_weighted_points(pdb, t, weighted_points); 
  
//   dsrpdb::Model model = pdb.model(0);
//   for (unsigned int i=0; i<model.number_of_chains(); i++) {
//     dsrpdb::Protein protein = model.chain(i);
//     for (typename dsrpdb::Protein::Atoms_iterator it= protein.atoms_begin(); 
// 	 it != protein.atoms_end();
// 	 ++it){
//       // NGHK: Set radius:
//       Weighted_point *wp = new Weighted_point(it->second.cartesian_coords(),1);
//       *weighted_points++ = *wp;
//    }
//   }
}
			    

#endif // EXTRACT_BALLS_FROM_PDB_H
