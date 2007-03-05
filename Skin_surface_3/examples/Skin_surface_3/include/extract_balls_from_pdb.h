// AUTHOR: Daniel Russel drussel@graphics.stanford.edu

// See http://graphics.stanford.edu/~drussel/pdb for documentation.
// 
// Feel free to contact me (Daniel, drussel@graphics.stanford.edu) 
// if you have any suggestions or patches. Thanks.

#ifndef EXTRACT_BALLS_FROM_PDB_H
#define EXTRACT_BALLS_FROM_PDB_H

#include <CGAL/PDB/PDB.h>
#include <CGAL/PDB/cgal.h>
#include <fstream>

template <class Traits, class OutputIterator>
void extract_balls_from_pdb(const char *filename, 
			    Traits const &t,
			    OutputIterator weighted_points) 
{
  std::ifstream in(filename);

  if (!in) {
    std::cerr << "Error opening input file " << filename << std::endl;
  }

  CGAL::PDB::PDB pdb(in, true);
  CGAL::PDB::all_weighted_points(pdb.model(0), t, weighted_points); 
}
			    

#endif // EXTRACT_BALLS_FROM_PDB_H
