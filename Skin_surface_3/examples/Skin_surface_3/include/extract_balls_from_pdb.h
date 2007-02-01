// AUTHOR: Daniel Russel drussel@graphics.stanford.edu

// See http://graphics.stanford.edu/~drussel/pdb for documentation.
// 
// Feel free to contact me (Daniel, drussel@graphics.stanford.edu) 
// if you have any suggestions or patches. Thanks.

#ifndef EXTRACT_BALLS_FROM_PDB_H
#define EXTRACT_BALLS_FROM_PDB_H

// include the source files, since the generated makefile 
// cannot link multiple object files
#include "dsrpdb/lib/pdb_utils.cc"
#include "dsrpdb/lib/Residue.cc"
#include "dsrpdb/lib/Protein.cc"
#include "dsrpdb/lib/Error_logger.cc"
#include "dsrpdb/lib/Protein_pdb.cc"
#include "dsrpdb/lib/Model.cc"
#include "dsrpdb/lib/PDB.cc"
#include "dsrpdb/lib/Residue_data.cc"
    
#include "dsrpdb/PDB.h"
#include "dsrpdb/geometry.h"
#include <fstream>

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
}
			    

#endif // EXTRACT_BALLS_FROM_PDB_H
