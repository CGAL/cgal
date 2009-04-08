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
#include <CGAL/PDB/range.h>
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
  CGAL::PDB::PDB p(in);
  //p.write(std::cout);
  std::ostringstream of;
	
  std::cout << "There are " << p.models().size() << " models." << std::endl;


  CGAL_PDB_FOREACH(const CGAL::PDB::PDB::Model_pair &m, p.models()) {
    std::cout << "Model " << m.key() << " has " 
	      << m.model().chains().size() << " chains" << std::endl;
	
    unsigned int na= CGAL::PDB::distance(m.model().atoms());
    std::cout << "There are " << na << " atoms" << std::endl;
    unsigned int nb= CGAL::PDB::distance(m.model().bonds());
    std::cout << "There are " << nb << " bonds" << std::endl;

    assert(na==1059);
    assert(nb==1062);

    unsigned int totaL_atoms=0;
    unsigned int total_bonds=0;
    CGAL_PDB_FOREACH(const CGAL::PDB::Model::Chain_pair &c, m.model().chains()) {
      std::cout << "Chain " << c.key() << " has " 
		<<c.chain().monomers().size() << " residues" 
		<< std::endl;
      totaL_atoms += CGAL::PDB::distance(c.chain().atoms());
      total_bonds += CGAL::PDB::distance(c.chain().bonds());
    }
    assert(CGAL::PDB::distance(m.model().atoms()) == totaL_atoms);
    assert(CGAL::PDB::distance(m.model().bonds()) == total_bonds);
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
