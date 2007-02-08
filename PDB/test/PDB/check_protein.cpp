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
along with the DSR PDB Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <CGAL/PDB/Protein.h>
#include <CGAL/PDB/geometry.h>
#include <CGAL/PDB/iterator.h>
#include <fstream>
#include <cassert>
#include <iterator>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include "check_equal.h"

int main(int argc, char *argv[]){
  //dsr::Residue res= dsr::Residue(dsr::Residue::VAL);
  //res.write(std::cout);
  
  std::ifstream in("data/check_protein.pdb");
  CGAL_PDB_NS::Protein p(in);
  //p.write(std::cout);
  std::ostringstream of;
  
  std::cout << "There are " << p.number_of_residues() << " residues." << std::endl;
  
  p.write_pdb(of);
  
  std::ifstream in2("data/check_protein.pdb");
  std::istringstream iss(of.str().c_str());
  check_equal(in2, iss);
  
  //p.dump(std::cout);

  // Demonstrate the geometry functions
  std::vector<CGAL_PDB_NS::Point> atms;
  CGAL_PDB_NS::coordinates(p.atoms_begin(), p.atoms_end(), std::back_inserter(atms));
  std::cout << "There are " << atms.size() << " atoms." << std::endl;
  assert(std::distance(p.atoms_begin(), p.atoms_end()) == atms.size());
  
  std::vector<CGAL_PDB_NS::Protein::Bond> bds(p.bonds_begin(), p.bonds_end());
  std::cout << "There are " << bds.size() << " bonds." << std::endl;
  

  std::vector<CGAL_PDB_NS::Point> bb;
  CGAL_PDB_NS::backbone_coordinates(p.atoms_begin(), p.atoms_end(), std::back_inserter(bb));
  assert(bb.size() < atms.size());
  {
    std::vector<std::pair<int,int> > bonds;
    std::vector<CGAL_PDB_NS::Point> points;
    CGAL_PDB_NS::coordinates_and_bonds(p, std::back_inserter(points), std::back_inserter(bonds));
    assert(bonds.size() == bds.size());
    assert(points.size() == atms.size());
  }
  {
    std::vector<CGAL_PDB_NS::Point> points;
    CGAL_PDB_NS::coordinates(p.atoms_begin(), p.atoms_end(), std::back_inserter(points));
    assert(points.size() == atms.size());
  }
  {
    std::vector<CGAL_PDB_NS::Point> points;
    std::copy(backbone_coordinates_begin(p),
	      backbone_coordinates_end(p),
	      std::back_inserter(points));
    assert(points.size() == bb.size());
  }
  {
    std::vector<std::pair<int,int> > bonds;
    std::vector<CGAL_PDB_NS::Point> points;
    CGAL_PDB_NS::simplified_coordinates_and_bonds(p, 
					     std::back_inserter(points), 
					     std::back_inserter(bonds));
    //std::cerr << bonds.size()  << " " << atms.size() << " " << points.size() << std::endl;
    assert(bonds.size() == points.size()-1);
    assert(points.size() <= 4*p.number_of_residues());
  }

  {
    for (CGAL_PDB_NS::Protein::Atoms_iterator it= p.atoms_begin(); it != p.atoms_end(); ++it){
      CGAL_PDB_NS::Atom::Index ind= it->second.index();
      CGAL_PDB_NS::Residue::Atom_label al= it->first;
      CGAL_PDB_NS::Residue& res= p.residue_containing_atom(ind);
      assert(&res >= &*p.residues_begin() 
	     && &res < &*p.residues_begin() + p.number_of_residues());
      assert(res.atom_label(ind) != CGAL_PDB_NS::Residue::AL_INVALID);
      assert(res.atom_label(ind) == al);
    }
  }

  {
    for (CGAL_PDB_NS::Protein::Bonds_iterator it= p.bonds_begin(); it != p.bonds_end(); ++it){
      CGAL_PDB_NS::Atom::Index ind0= it->first;
      CGAL_PDB_NS::Atom::Index ind1= it->second;
      assert(p.atom(ind0) != p.atom(CGAL_PDB_NS::Atom::Index()));
      assert(p.atom(ind1) != p.atom(CGAL_PDB_NS::Atom::Index()));
      int ir0= p.residue_containing_atom(ind0).index().to_index();
      int ir1= p.residue_containing_atom(ind1).index().to_index();
      int diff = ir1 - ir0;
      assert(diff==0 || diff ==1);
    }
  }

  return return_code__;
}
