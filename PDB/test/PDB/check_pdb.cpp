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
#include <sstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include "check_equal.h"


int main(int , char *[]){
  //dsr::Residueres= dsr::Residue(dsr::Residue::VAL);
  //res.write(std::cout);
  std::string	argv1="data/check_pdb.pdb";
  std::ifstream in(argv1.c_str());
  CGAL_PDB_NS::PDB p(in);
  //p.write(std::cout);
  std::ostringstream of;
	
  std::cout << "There are " << p.number_of_models() << " models." << std::endl;


  for (CGAL_PDB_NS::PDB::Model_iterator it= p.models_begin();
       it != p.models_end(); ++it){
    const CGAL_PDB_NS::Model &m= it->model();
    std::cout << "Model " << it->key() << " has " 
	      << m.number_of_chains() << " chains" << std::endl;
	
    unsigned int na= std::distance(m.atoms_begin(), m.atoms_end());
    std::cout << "There are " << na << " atoms" << std::endl;
    unsigned int nb= std::distance(m.bonds_begin(), m.bonds_end());
    std::cout << "There are " << nb << " bonds" << std::endl;

    assert(na==1059);
    assert(nb==1062);

    unsigned int totaL_atoms=0;
    unsigned int total_bonds=0;
    for (CGAL_PDB_NS::Model::Chain_const_iterator it= m.chains_begin(); it != m.chains_end(); ++it){
      std::cout << "Chain " << it->key() << " has " 
		<<it->chain().number_of_monomers() << " residues" 
		<< std::endl;
      totaL_atoms += std::distance(it->chain().atoms_begin(), it->chain().atoms_end());
      total_bonds += std::distance(it->chain().bonds_begin(), it->chain().bonds_end());
    }
    assert(std::distance(m.atoms_begin(), m.atoms_end()) == totaL_atoms);
    assert(std::distance(m.bonds_begin(), m.bonds_end()) == total_bonds);
  }
	
  p.write(of);
  std::ofstream tmpfile("/tmp/out.pdb");
  p.write(tmpfile);
  //of << std::fflush;
	
  std::ifstream nif(argv1.c_str());
  std::istringstream ofis(of.str().c_str());
  check_equal(nif, ofis);
	
  return return_code__;
}
