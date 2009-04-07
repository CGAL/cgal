/* Copyright 2004
Stanford University

This file is part of the DSR PDB Library.

The DSR PDB Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The DSR PDB Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the DSR PDB Library; see the file LICENSE.LGPL.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA. */
#include <CGAL/PDB/PDB.h>
#include <fstream>
#include <cassert>

#include "include/check_equal.h"


int main(int , char *[]){
  //dsr::Residue res= dsr::Residue(dsr::Residue::VAL);
  //res.write(std::cout);
  //assert(argc==3);
  CGAL_SET_LOG_LEVEL(CGAL::Log::LOTS);
  std::ifstream in("data/check_atom_hetatom.pdb");
  CGAL::PDB::PDB p(in);
  //p.write(std::cout);
  std::ostringstream of;
  
  std::cout << "There are " << p.models().size() << " models." << std::endl;
  
  CGAL_PDB_FOREACH(const CGAL::PDB::PDB::Model_pair& m,
		           p.models()) {
    std::cout << "Model " << m.key() << " has " << m.model().chains().size() 
	      << " chains" << " and " << m.model().heterogens().size() 
	      << " heterogens" << std::endl;
    CGAL_PDB_FOREACH(const CGAL::PDB::Model::Chain_pair &c,
					m.model().chains()) {
	std::cout << "Chain " << c.key() << " has " 
		<< c.chain().monomers().size() << " residues" << std::endl;
    }
  }
  p.write(of);

  return EXIT_SUCCESS;
}
