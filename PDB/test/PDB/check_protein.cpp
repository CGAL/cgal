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

#include <CGAL/PDB/Chain.h>
#include <CGAL/PDB/iterator.h>
#include <fstream>
#include <CGAL/Testsuite/assert.h>
#include <iterator>
#include <CGAL/Tools/Log.h>

#include "check_equal.h"

int main(int , char *[]){
  //dsr::Residue res= dsr::Residue(dsr::Residue::VAL);
  //res.write(std::cout);
  
  std::ifstream in("data/check_protein.pdb");
  CGAL_PDB_NS::Chain p(in);
  //p.write(std::cout);
  std::ostringstream of;
  
  std::cout << "There are " << p.number_of_monomers() << " residues." << std::endl;
  
  p.write_pdb(of);
  
  /*std::ofstream out("/tmp/out.pdb");
    p.write_pdb(out);*/
 
    std::ifstream in2("data/check_protein.pdb");
  std::istringstream iss(of.str().c_str());
  check_equal(in2, iss);
  
  unsigned int na= std::distance(p.atoms_begin(), p.atoms_end());
  std::cout << "There are " << na << " atoms" << std::endl;
  unsigned int nb= std::distance(p.bonds_begin(), p.bonds_end());
  std::cout << "There are " << nb << " bonds" << std::endl;
  CGAL_assert_equal(na, 1221);
  CGAL_assert_equal(nb, 1248);

  return return_code__;
}
